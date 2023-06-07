# -*- coding: utf-8 -*-
"""
Organism object
It allocates the full data structure
"""
import math
import random
import time
import numpy as np
import copy
from .placement_object import PlacementObject
from .shape_object import ShapeObject
from .connector_object import norm_cdf
from .connector_object import g
import _multiplacement

class OrganismObject:
    """Organism object
       The Organism object essentially contains two vectors:
       - A vector of connector objects
       - A vector of recognizer objects
       The order of the elements in these vectors determine, implicitly,
       the connections between the elements (i.e. what recognizers a
       particular connector is connecting)
    """

    def __init__(self, _id: int, conf: dict, max_pssm_length: int) -> None:
        """Organism constructor

        Args:
            _id: Organism identifier (assigned by factory)
            conf: Organism-specific configuration values from JSON file
            self.pwm_length: Maximum column size of the pssm recognizer
	    [will prevent mutations from going over max pssm length]
        """
        self._id = _id
        
        # Instantiate recognizer and connector vectors
        self.recognizers = []
        self.connectors = []
        self.recognizers_flat = []
        self.connectors_scores_flat = []
        self.recognizer_lengths = []
        self.recognizer_types = ""
        self.recognizer_models = []
        self.recognizer_bin_edges = []
        self.recognizer_bin_nums = []
        self.minimum_length = 0
        # assign organism-specific parameters
	
        # whether fitness is computed over sequences as sum or average
        self.is_precomputed = conf["PRECOMPUTE"]
        self.cumulative_fit_method = conf["CUMULATIVE_FIT_METHOD"]
        
        # energy thresholding parameters (to prevent negative energy values)
        self.energy_threshold_method = conf["ENERGY_THRESHOLD_METHOD"]
        self.energy_threshold_value = conf["ENERGY_THRESHOLD_PARAM"]
        
        # probability of replacing PSSM by random PSSM
        self.mutate_probability_substitute_pssm = conf[
            "MUTATE_PROBABILITY_SUBSTITUTE_PSSM"
        ]
        
        # type of indel operator: 
        # - blind (preserves connectors' mus and sigmas)
        # - intelligent (preserves relative distances between recognizers) 
        self.insertion_method = conf[
            "INSERTION_METHOD"
        ]
        self.deletion_method = conf[
            "DELETION_METHOD"
        ]
	
        #probability of deleting/inserting a recognizer
        self.mutate_probability_delete_recognizer = conf[
            "MUTATE_PROBABILITY_DELETE_RECOGNIZER"
        ]
        self.mutate_probability_insert_recognizer = conf[
            "MUTATE_PROBABILITY_INSERT_RECOGNIZER"
        ]
	
        # probability of mutating a node
        self.mutate_probability_node_mutation = conf[
            "MUTATE_PROBABILITY_NODE_MUTATION"
        ]
        
        # min and maximum number of nodes allowed
        self.min_nodes = conf["MIN_NODES"]
        self.max_nodes = conf["MAX_NODES"]
        
        # maximum length of PSSMs allowed
        self.max_pssm_length = max_pssm_length
        
        # Map used by the placement algorithm
        # The list maps each row of the matrix of the placement scores onto a
        # column of a PSSM: each row is assigned a [pssm_idx, column_idx]
        self.row_to_pssm = []
        
        # Dictionary storing information about how the organism has to be
        # assembled by the recombination process. All the values are
        # initialized as None.
        # The "recognizers" key maps to the list of required recognizers, from
        # left to right. Same thing for the "connectors" key. The "p1" and "p2"
        # keys map to the ID number of parent 1 and parent 2, respectively.
        self.assembly_instructions = {'p1': None,  # ID number of the first parent
                                      'p2': None,  # ID number of the second parent
                                      'recognizers': None,
                                      'connectors': None}
    
    def set_assembly_instructions(self, aligned_repres, connectors_table, p1_id, p2_id):
        '''
        Sets the self.assembly_instructions attribute.
        '''
        self.assembly_instructions['p1'] = p1_id
        self.assembly_instructions['p2'] = p2_id
        self.set_recogs_assembly_instructions(aligned_repres)
        self.set_connectors_assembly_instructions(aligned_repres, connectors_table)
    
    def set_recogs_assembly_instructions(self, aligned_repres):
        '''
        Compiles the list of recognizers in the self.assembly_instructions
        attribute.
        '''
        child_recog_names = [item for item in aligned_repres if item != '-']
        recogs_list = []
        for item in child_recog_names:
            if item not in recogs_list:
                recogs_list.append(item)
        self.assembly_instructions['recognizers'] = recogs_list
    
    def set_connectors_assembly_instructions(self, aligned_repres, connectors_table):
        '''
        Compiles the list of connectors in the self.assembly_instructions
        attribute.
        '''
        required_connectors = self.annotate_required_connectors(aligned_repres)
        connectors_list = []
        for (left, right) in required_connectors:
            # What suitable connectors are available to cover that span
            available_connectors = connectors_table[left][right]
            
            # If there's already a suitable connector
            if len(available_connectors) > 0:
                # There can be one or two (from the two parents) equivalent
                # connectors available. Randomly pick one.
                connector_name = available_connectors.pop(
                    random.randint(0, len(available_connectors)-1)
                )
            
            # If no connector is available to cover that span, make an appropriate one
            else:
                connector_name = 'synth_' + str(left) + '_' + str(right)
            
            connectors_list.append(connector_name)
        self.assembly_instructions['connectors'] = connectors_list
    
    def get_parent1_parent2_ratio(self):
        '''
        This function is used on newly assembled children, obtained by
        recombination between two parent organisms. The assembly_instructions
        attribute reveals how many nodes came from each parent. This method
        counts how many nodes come from parent1, and how many from parent2.
        Then it returns the ratio between these two counts.
        '''
        
        from_p1 = 0  # How many nodes from parent1
        from_p2 = 0  # How many nodes from parent2
        
        # Count recognizers from each parent
        for node_name in self.assembly_instructions['recognizers']:
            if node_name[:3] == 'p1_':
                from_p1 += 1
            elif node_name[:3] == 'p2_':
                from_p2 += 1
                
        # Count connectors from each parent
        for node_name in self.assembly_instructions['connectors']:
            if node_name[:3] == 'p1_':
                from_p1 += 1
            elif node_name[:3] == 'p2_':
                from_p2 += 1
        
        # Compute parent1 / parent2 ratio
        if from_p2 == 0:  # Avoid 0-division error
            p1_p2_ratio = np.inf
        else:
            p1_p2_ratio = from_p1 / from_p2
        
        return p1_p2_ratio
    
    def annotate_required_connectors(self, child_repres):
        '''
        Returns the list of the connectors' spans required to join the
        recognizers, given the representation of the organism (a newly
        generated child). Each element of the list is a tuple of two element,
        where the first one is the index position of the recognizer to the left,
        while the secind one is the index position of the recognizer to the right.
        
        EXAMPLE:
        In this representation
        
            p1_0    p2_0    -       p1_2
        
        A connector is needed to link 'p1_0', which is at index position 0, to
        'p2_0', which is at index position 1. Another connector is needed to
        link 'p2_0' with 'p1_2', which is at index position 3. Therefore, the
        returned list of required connectors would be
        
            [ (0, 1), (1, 3) ]
        
        In case of repetitions (due to a recognizer overlapping with more than
        one recog in the alignment of the organisms), no recognizer is needed.
        
        EXAMPLE:
        In this representation
        
            p1_0    -       p1_1    p1_1    p2_3
        
        The returned list of required connectors would be
        
            [ (0, 2), (2, 3) ]
        
        No connector is needed for index 2 to index 3: positions that are both
        occupied by the same recognizer ('p1_1').
        
        '''
        # Indexes where there is a recognizer
        recogs_indexes = []
        for idx in range(len(child_repres)):
            if child_repres[idx] != '-':
                recogs_indexes.append(idx)
        
        required_connectors = []
        for i in range(len(recogs_indexes)-1):
            
            left_recog_index = recogs_indexes[i]
            right_recog_index = recogs_indexes[i+1]
            
            left_recog_name = child_repres[left_recog_index]
            right_recog_name = child_repres[right_recog_index]
            
            # No connector is needed to link a recognizer with itself (see
            # function's docstring).
            if left_recog_name != right_recog_name:
                connector_span = (left_recog_index, right_recog_index)
                required_connectors.append(connector_span)
        
        return required_connectors
    
    def set_row_to_pssm(self):
        """row_to_pssm is an attribute that maps each row of the alignment
           matrix to the index (within org) of the pssm recognizer, and the
           column of the pssm that applies to that row.
           In the alignment matrix, the rows correspond all consecutively to 
           pssm positions (i.e. columns).
           This attribute allows to go from row directly to a pair of indices
           denoting the pssm number and the position within the pssm that maps
           to it.
           
           Call by factory upon generation of organism, and also after any
           mutations that might change the size of the pssm, or the column order,
           or their number.
        """
        
        pssm_list = self.recognizers
        
        # Initialize list
        # (first row in the placement matrix doesn't represent any PSSM position)
        row_to_pssm_list = [[None, None]]
        
        # Fill the list
        for i in range(len(pssm_list)):
            for j in range(pssm_list[i].length):
                row_to_pssm_list.append([i, j])
        
        # Token row
        # (this ensures that the last row will be considered the row of the last
        # position of a PSSM when calling self.is_last() on it)
        row_to_pssm_list.append([None, 0])
        
        self.row_to_pssm = row_to_pssm_list

    def get_id(self) -> int:
        """Getter _id

        Returns:
            _id of the organism
        """
        return self._id

    def set_id(self, _id: int) -> None:
        """Setter _id

        Args:
            _id: ID to to set in the organism
        """
        self._id = _id
    
    def compress_connectors(self, total_span_avg, total_span_var):
        '''
        When a new recognizer is inserted, two connectors do the work that was
        previously done by one connector. That connector is called the
        "old connector". Both the old connector and the newly inserted connector
        are compressed (mu and sigma are shrinked) so that together they model an
        overall distance equivalently to how the old connector was modeling it
        (assuming they are two independent normally distributed random variables).
        
        Args:
            total_span_avg:
                the sum of the two mu values should be equal to `total_span_avg`
            total_span_var:
                the sum of the two variances should be equal to `total_span_var`
        
        Returns:
            the adjusted parameters (mu and sigma) for the two connectors
        '''
        
        # Randomly chose two mu values such that the sum is old_mu
        mu_left = random.random() * total_span_avg
        mu_right = total_span_avg - mu_left
        # Randomly chose two var values such that the sum is old_var
        var_left = random.random() * total_span_var
        var_right = total_span_var - var_left
        # sigma values
        sigma_left = np.sqrt(var_left)
        sigma_right = np.sqrt(var_right)
        return [(mu_left, sigma_left), (mu_right, sigma_right)]
    
    def mutate(self, org_factory) -> None:
        """Mutates an organism based on JSON configured probabilities

        Args:
            org_factory (OrganismFactory): Factory of organisms and node
                                           components
        """

        
        # Deletion
        if random.random() < self.mutate_probability_delete_recognizer:
            # Delete a recognizer (and a connector)
            self.delete_a_recognizer()
        
        # Insertion
        if random.random() < self.mutate_probability_insert_recognizer:
            # Insert a recognizer (and a connector)
            self.insert_a_recognizer(org_factory)
        
        # Node-mutations
        self.mutate_nodes(org_factory)
        
        # After applying mutations, order/columns/number of pssm's may have
        # changed, so we call the set_row_to_pssm to set their mapping on the
        # alignment matrix anew
        self.set_row_to_pssm()
        self.flatten()
    
    def delete_a_recognizer(self):
        '''
        A randomly selected recognizer is deleted.
        
        Connectors are adjusted as needed:
            
            - If the deleted recognizer was a terminal node, its adjacent
              connector is also deleted.
              
            - If the deleted recognizer was an internal node, its two adjacent
              connectors are replaced by a single connector.
              
              Note:
              Depending on the `deletion_method`, specified in the config, this
              replacement can occur by randomly chosing which connector to keep
              ("blind" method) or by merging their normal random variables into
              a single new random variable with a mu and sigma that model the
              spacer equivalently to the sum of the two pre-existing connectors
              ("intelligent" method).
        '''
        n_recognizers = self.count_recognizers()
        # if this is a single-node organism, skip deletion
        if n_recognizers > 1:
            
        	# "blind" method: the remaining connector is left unchanged
            if self.deletion_method == "blind":
                # Choose randomly the recognizer to be deleted
                recognizer_idx = random.randint(0, n_recognizers - 1)
                if recognizer_idx == 0:
                    # If the recognizer to delete is the first, the connector to
                    # delete is the one to the right
                    # ( same index: connector_idx = recognizer_idx = 0 )
                    connector_idx = recognizer_idx
                elif recognizer_idx == n_recognizers - 1:
                    # If the recognizer to delete is the last, the connector to
                    # delete is the one to the left
                    # ( different index: connector_idx = recognizer_idx - 1 )
                    connector_idx = recognizer_idx - 1
                else:
                    # if the recognizer to delete is in not a terminal recognizer
                    # of the chain, the parent connector to be deleted (left/riht)
                    # is chosen randomly
                    if random.random() < 0.5:
                        connector_idx = recognizer_idx
                    else:
                        connector_idx = recognizer_idx - 1
                
				# recreate arrays of recognizers and connectors
				# skipping the recgonizer/connector selected for deletion
                new_recognizers = (self.recognizers[:recognizer_idx] +
                                   self.recognizers[recognizer_idx + 1:])
                new_connectors = (self.connectors[:connector_idx] +
                                   self.connectors[connector_idx + 1:])
        		# assign new arrays
                self.recognizers = new_recognizers
                self.connectors = new_connectors
            
 			# "intelligent" method: a new connector is created by merging
 			# the connectors on the sides of the deleted recognizer into a
            # single larger connector
            elif self.deletion_method == "intelligent":
                # Choose randomly the recognizer to be deleted
                recognizer_idx = random.randint(0, n_recognizers - 1)
                if recognizer_idx == 0:
                    # If the recognizer to delete is the first, the connector to
                    # delete is the one to the right
                    # ( same index: connector_idx = recognizer_idx = 0 )
                    connector_idx = recognizer_idx
                elif recognizer_idx == n_recognizers - 1:
                    # If the recognizer to delete is the last, the connector to
                    # delete is the one to the left
                    # ( different index: connector_idx = recognizer_idx - 1 )
                    connector_idx = recognizer_idx - 1
                else:
                    # if the recognizer to delete is in not a terminal recognizer
                    # of the chain, the parent connector to be deleted (left/right)
                    # is chosen randomly
                    if random.random() < 0.5:
                        # Index of the parent connector that will be deleted
                        connector_idx = recognizer_idx
                        # Index of the parent connector that will be adjusted
                        connector_to_stretch = recognizer_idx - 1
                    else:
                        # Index of the parent connector that will be deleted
                        connector_idx = recognizer_idx - 1
                        # Index of the parent connector that will be adjusted
                        connector_to_stretch = recognizer_idx
                    
                    # Adjust parameters of the neighbour connector
                    ''' The parent connector that is not deleted is modified,
                    so that it can span the gap left by the deletion, without
                    heavily affecting the placement of the nodes to the sides of
                    the deletion point.'''
                    
                    # Adjust MU
                    '''Thus, we need to add to the mu of the remaining connector
                    also the mu of the deleted connector and the legth of the
                    deleted recognizer.'''
                    adj_mu = (self.connectors[connector_to_stretch]._mu +
                              self.connectors[connector_idx]._mu +
                              self.recognizers[recognizer_idx].length)
                    
                    # Adjust SIGMA
                    ''' Variances are additive (under assumption of independence
                    between the two random variables). Therefore to adjust the
                    variance we need to add the variance of the deleted
                    connector (the length of the deleted recognizer has no
                    variance). Therefore, its standard deviation (sigma) becomes
                    the square root of the sum of the squares of the sigmas.'''
                    adj_sigma = (self.connectors[connector_to_stretch]._sigma ** 2 +
                                self.connectors[connector_idx]._sigma ** 2)**(1/2)
                    # set new mu and new sigma
                    self.connectors[connector_to_stretch].set_mu(adj_mu)
                    self.connectors[connector_to_stretch].set_sigma(adj_sigma)
                
				# recreate arrays of recognizers and connectors
				# skipping the recgonizer/connector selected for deletion					
                new_recognizers = (self.recognizers[:recognizer_idx] +
                                   self.recognizers[recognizer_idx + 1:])
                new_connectors = (self.connectors[:connector_idx] +
                                   self.connectors[connector_idx + 1:])
                
        		# assign new arrays
                self.recognizers = new_recognizers
                self.connectors = new_connectors
            
            else:
                raise ValueError('DELETION_METHOD in the config should be ' +
                                 '"blind" or "intelligent".')
    
    def insert_a_recognizer(self, org_factory):
        '''
        A new recognizer is inserted in a random position of the chain.
        
        Connectors are adjusted as needed:
            
            - If the inserted recognizer is a terminal node, a new random
              connector is used to link it to the chain.
              
            - If the recognizer is inserted internally in the chain, the
              connector modeling the spacer where the insertion happened is
              replaced by two connectors.
              
              Note:
              Depending on the `insertion_method`, specified in the config, this
              replacement can occur by adding a random connector to the
              pre-existing one ("blind" method) or by redefining the parameters
              of the two connectors (the pre-existing one and the newly generated
              one) so that the sum of the two normal random variables is equivalent
              to the normal random variable originally specified by the pre-existing
              connector before the insertion ("intelligent" method).
        '''
        
        # instantiate the new recognizer and the new connector
        new_connector = org_factory.create_connector()
        new_recognizer = org_factory.create_recognizer()
       
        # "blind" method: pre-existing connector linked to the new inserted
        # recognizer is left unmodified
        if self.insertion_method == "blind":
            
            n_recognizers = self.count_recognizers()
            # Choose randomly the recognizer next to which the insertion
            # is going to occur
            recognizer_idx = random.randint(0, n_recognizers - 1)
            # Choose randomly whether the insertion is going to be to the
            # left or to the right of the considered recognizer
            
            if random.random() < 0.5:  # Insertion occurs to the left
                # First connector after insertion point and first
                # recognizer after insertion point
                connector_idx = recognizer_idx
                
            else:  # Insertion occurs to the right
                # First connector after insertion point
                connector_idx = recognizer_idx
                # First recognizer after insertion point
                recognizer_idx += 1
                
            # recreate arrays of connectors/recognizers, adding
            # the newly minted recognizer+connector
            new_recognizers = (self.recognizers[:recognizer_idx] +
                               [new_recognizer] +
                               self.recognizers[recognizer_idx:])
            new_connectors = (self.connectors[:connector_idx] +
                               [new_connector] +
                               self.connectors[connector_idx:])
            
            # assign new connector/recognizer arrays
            self.recognizers = new_recognizers
            self.connectors = new_connectors
            
            
        # "intelligent" method: the pre-existing connector and the new one are
        # "shrinked" so that together they model a spacer "equivalent" to
        # the one modeled by the connector previously occupying the space
        # where the new recognizer has been inserted
        elif self.insertion_method == "intelligent":
            
            n_recognizers = self.count_recognizers()
            # Choose randomly the recognizer next to which the insertion
            # is going to occur
            recognizer_idx = random.randint(0, n_recognizers - 1)
            
            # set no compression as default (for terminal insertion cases)
            connector_to_compress = None
            
            # Choose randomly whether the insertion is going to be to the
            # left or to the right of the considered recognizer
            if random.random() < 0.5:  # Insertion occurs to the left
                # First connector after insertion point and first
                # recognizer after insertion point
                connector_idx = recognizer_idx
                
            	# if the new recognizer is NOT the first in the chain
                if recognizer_idx != 0:
                    connector_to_compress = recognizer_idx - 1
                    # (No compression is required if insertion occurs to
                    # the left of the first recognizer of the chain)
            
            else:  # Insertion occurs to the right
                # First connector after insertion point
                connector_idx = recognizer_idx
                
            	# if the new recognizer is NOT the last in the chain
                if recognizer_idx != n_recognizers - 1:
                    connector_to_compress = recognizer_idx
                    # (No compression is required if insertion occurs to
                    # the right of the last recognizer of the chain)
                    
                    # First recognizer after insertion point
                    recognizer_idx += 1
            
            # if connector needs to be "compressed" (not a terminal insertion)
            if connector_to_compress != None:
                
                '''Ideally, we would like the sum of the mus of the two
                connectors (the old one and the inserted one) + the length
                of the inserted recognizer to be equal to the gap spanned
                by the pre-existing connector, so that the insertion
                doesn't heavily affect the placement of the nodes to the
                sides of the insertion point. So we need to shrink the
                mus and sigmas accordingly.'''
                
                # Adjust MUs approach
                ''' If the legth of the inserted recognizer alone is larger
                than the gap (i.e. larger than the mu of the already present
                connector) we can't reach the goal entirely and the best
                thing we can do is to keep the mus of the two connectors
                maximally shrinked. Otherwise, the mus of the connectors
                will be scaled so that their sum is equal to the gap (mu of
                the old connector) minus the length of the inserted
                recognizer.'''
                
                # Adjust SIGMAs approach
                ''' Variances are additive (under assumption of independence
                between the two random variables). Therefore, the overall
                variance will be the sum of the variances of the two
                connectors. The variances will be scaled so that their sum
                is equal to the variance of the connector that was already
                present before the insertion. The standard deviations,
                which are the square root of the variances, will be adjusted
                accordingly.'''
                
                # Define what the sum of the two mu values should be (expected_sum_mus)
                mu_old_connector = self.connectors[connector_to_compress]._mu
                expected_sum_mus = mu_old_connector - new_recognizer.length
                # If the inserted recognizer alone is larger than the gap
                if expected_sum_mus < 0:
                    # best thing we can do is to maximally shrink mus
                    expected_sum_mus = 0
                # Define what the variance of the total span should be (it
                # should be equal to the variance of the pre-existing connector)
                expected_var = self.connectors[connector_to_compress]._sigma ** 2
                # Updated values for the two connectors that need compression
                preexisting, inserted = self.compress_connectors(expected_sum_mus, expected_var)
                # Update parameters of pre-existing connector
                self.connectors[connector_to_compress].set_mu(preexisting[0])
                self.connectors[connector_to_compress].set_sigma(preexisting[1])
                # Update paramters of inserted connector
                new_connector.set_mu(inserted[0])
                new_connector.set_sigma(inserted[1])
                    
            # recreate arrays of connectors/recognizers, adding
            # the newly minted recognizer+connector and containing
            # also the "compressed" existing connector
            new_recognizers = (self.recognizers[:recognizer_idx] +
                               [new_recognizer] +
                               self.recognizers[recognizer_idx:])
            new_connectors = (self.connectors[:connector_idx] +
                               [new_connector] +
                               self.connectors[connector_idx:])
        
        else:
            raise ValueError('INSERTION_METHOD in the config should be ' +
                             '"blind" or "intelligent".')
        
        # assign new connector/recognizer arrays
        self.recognizers = new_recognizers
        self.connectors = new_connectors
    
    def mutate_nodes(self, org_factory):
        '''
        Triggers node-specific mutations on the nodes of the organism.
        
        If `mutate_probability_node_mutation` is specified, only one node will
        become prone to mutations, with a probability equal to
        `mutate_probability_node_mutation`. Otherwise, all the nodes are always
        prone to mutations.
        
        When a node is prone to mutations it doesn't mean it will necessarily
        be mutated. Each type of node mutation can occur with its own specific
        proability (specified in the config).
        '''
        
        # If MUTATE_PROBABILITY_NODE_MUTATION is set to real value, a single
        # node is selected. If set to null, all nodes are prone to mutations
        if self.mutate_probability_node_mutation:
            # Mutate a random node
            if random.random() < self.mutate_probability_node_mutation:
                
                n_nodes = self.count_nodes()
                random_node_idx = random.randint(0, n_nodes - 1)
                if random_node_idx < self.count_recognizers():
                    # mutate a recognizer
                    moved_pssm_bounds = self.recognizers[random_node_idx].mutate(org_factory)
                    # Adjust adjacent gaps if needed
                    if moved_pssm_bounds:
                        self.adjust_gaps_after_pssm_bounds_displacement(
                            random_node_idx, moved_pssm_bounds)
                else:
                    # mutate a connector
                    connector_idx = random_node_idx - self.count_recognizers()
                    self.connectors[connector_idx].mutate(org_factory)
        
        else:  # All nodes are prone to mutations
            
            # Recognizers
            for i in range(self.count_recognizers()):
                # Mutate recognizer with index i
                moved_pssm_bounds = self.recognizers[i].mutate(org_factory)
                
                # Adjust adjacent gaps if needed
                if moved_pssm_bounds:
                    self.adjust_gaps_after_pssm_bounds_displacement(i, moved_pssm_bounds)
            
            # Connectors
            for connector in self.connectors:
                connector.mutate(org_factory)

    def adjust_gaps_after_pssm_bounds_displacement(self, pssm_index,
                                                   pssm_displacement_code):
        '''
        For well-adapted organisms (organisms that have a consistent placement
        strategy on the positive set) some PSSM mutations will produce a slight
        change in the DNA positions occupied by the PSSMs. If the PSSM gets a
        new column added to the left, and if there's a gap to the left of that
        PSSM, the connector modeling that gap should decrease its mu by 1, so
        that it mantains its optimality. Otherwise, the benefit of a good extra
        PSSM column could be overcomed by the downsides of having the average
        gap size off by 1 unit in the adjacent connector.
        
        The PSSM mutations to be taken into account for this problem are:
            - increase/decrease PSSM size
            - right/left shift of the PSSM
        Those are the mutation that can produce some displacement of the PSSM's
        bounds (the left and the right bounds). When they're called, a code is
        generated, that keeps track of how the bounds have changed. This code
        can be used as the second parameter of this function. The purpose of
        this function is to adjust connectors according to those codes, so that
        PSSMs are free to change size or shift without damaging adjacent
        connectors.

        Parameters
        ----------
        pssm_index : int
            Index of the PSSM that has gone through some bound displacement.
        pssm_displacement_code : list
            It's a list of two integers. The first integer is the shift
            occurred the left bound of the PSSM. The second integer is the
            shift occurred to its right bound.
            A positive shift means that the bound moved to the right, while a
            negative shift means a shift to the left.
            These codes are generated when mutating a PSSM with its mutate()
            method.

        '''
        
        left_bound_shift, right_bound_shift = pssm_displacement_code
        
        # If the LEFT bound has changed
        if left_bound_shift != 0:
            # Adjust connector to the left
            conn_index = pssm_index - 1
            
            # There's no connector the left of a PSSM that's the first
            # recognizer of the organism (no connector with index -1), so check
            # that the index is not -1.
            if conn_index >= 0:
                # Adjust mu
                new_mu = self.connectors[conn_index]._mu + left_bound_shift
                self.connectors[conn_index].set_mu(max(0, new_mu))
                
                # Update connector's PDF and CDF values
                self.connectors[conn_index].set_precomputed_pdfs_cdfs()
        
        # If the RIGHT bound has changed
        if right_bound_shift != 0:
            # Adjust connector to the right
            conn_index = pssm_index
            
            # There's no connector the right of a PSSM that's the last
            # recognizer of the organism, so check that the index is not equal
            # to the number of connectors
            if conn_index < self.count_connectors():
                # Adjust mu
                new_mu = self.connectors[conn_index]._mu - right_bound_shift
                self.connectors[conn_index].set_mu(max(0, new_mu))
                
                # Update connector's PDF and CDF values
                self.connectors[conn_index].set_precomputed_pdfs_cdfs()
        self.flatten()

    def get_M_and_SE_on_DNA_set(self, dna_set, method, gamma):
        '''
        Returns the mean and standard error to be used as input for the
        statistical test that we use as fitness function. Several variants of
        this test are possible. Depending on the variant of choice, the mean
        and standard error could be computed on a trimmed set.
        
        A lower bound of 1 unit is applied to the standard deviation (used to
        calculate the standard error). Being more consistent than that on a set
        of DNA sequences doesn't further improve fitness. This means that if the
        difference between averages (P - N) is small, the organism can not
        improve fitness indefinitely by only reducing the variances. [this would
        also make sense under the assumptino that discrimination happens in a
        noise-prone system where a very small P-N would be unreliable and
        ultimately wouldn't work].
        
        Methods description
        -------------------
        
        Welch:
            NO TRIMMING
        Yuen:
            TRIMMED MEAN.
            gamma*100 is the percentage of discarded elements from each tail of
            the distribution (both left and right)
        Trim-left-M:
            TRIMMED MEAN.
            gamma*100 is the percentage of discarded elements from the left tail of
            the distribution
        Trim-left-M-SD:
            TRIMMED MEAN AND TRIMMED STANDARD ERROR.
            gamma*100 is the percentage of discarded elements from the left tail of
            the distribution
        
        Parameters
        ----------
        dna_set : list
            Array of DNA sequences.
        method : str
            The fitness function of choice.
        gamma : float
            Number between 0 and 0.5.
        '''
        energy_scores = self.get_binding_energies(dna_set)
        n_to_trim = round(gamma * len(energy_scores))  # defalut gamma: 0.2
        # Make sure that `n_to_trim` doesn't exceed the maximum due to rounding
        # (when gamma=0.5 and len(energy_scores) is an odd number)
        # We should never trim more than half of the elements from one side.
        if n_to_trim > len(energy_scores)/2:
            n_to_trim = int(len(energy_scores)/2)
        
        if n_to_trim == 0:
            mean = np.mean(energy_scores)
            stdev = np.std(energy_scores)
            stdev = max(1, stdev)  # Lower bound to sigma (see docstring)
            sterr = stdev / np.sqrt(len(energy_scores))
        if method == "Welch":
            mean = np.mean(energy_scores)
            stdev = np.std(energy_scores)
            stdev = max(1, stdev)  # Lower bound to sigma (see docstring)
            sterr = stdev / np.sqrt(len(energy_scores))
        elif method == "Yuen":
            energy_scores.sort()
            energy_scores_trimmed = energy_scores[n_to_trim:-n_to_trim]
            mean = np.mean(energy_scores_trimmed)
            stdev = np.std(energy_scores)
            stdev = max(1, stdev)  # Lower bound to sigma (see docstring)
            sterr = stdev / np.sqrt(len(energy_scores))
        elif method == "Trim-left-M":
            energy_scores.sort()
            energy_scores_trimmed = energy_scores[n_to_trim:]
            mean = np.mean(energy_scores_trimmed)
            stdev = np.std(energy_scores)
            stdev = max(1, stdev)  # Lower bound to sigma (see docstring)
            sterr = stdev / np.sqrt(len(energy_scores))
        elif method == "Trim-left-M-SD":
            energy_scores.sort()
            energy_scores_trimmed = energy_scores[n_to_trim:]
            mean = np.mean(energy_scores_trimmed)
            stdev = np.std(energy_scores_trimmed)
            stdev = max(1, stdev)  # Lower bound to sigma (see docstring)
            sterr = stdev / np.sqrt(len(energy_scores_trimmed))
        else:
            raise ValueError('Unknown method "' + method + '". Please check ' +
                             'the FITNESS_FUNCTION paramter in the config file.')
        return mean, sterr
    
    def get_fitness(self, pos_set, neg_set, method, gamma=0.2):
        '''
        Returns the fitness of the organism.

        Parameters
        ----------
        pos_set : list
            The positive set of DNA sequences.
        neg_set : list
            The negative set of DNA sequences.
        method : str
            The fitness function of choice.
        gamma : float, optional
            Number between 0 and 0.5, used for trimming (used if required by
            the fitness function). The default is 0.2.
        '''
        M_p, SE_p = self.get_M_and_SE_on_DNA_set(pos_set, method, gamma)
        M_n, SE_n = self.get_M_and_SE_on_DNA_set(neg_set, method, gamma)
        return (M_p - M_n) / (SE_p**2 + SE_n**2)**(1/2)
    
    def get_binding_energies(self, dna_set):
        '''
        Returns a list of binding energies for all the sequences in `dna_set`.
        '''
        binding_energies = []
        for seq in dna_set:
            placement = self.get_placement(seq)
            binding_energies.append(placement.energy)
        return binding_energies

    def count_nodes(self) -> int:
        '''
        Returns the number of nodes of the organism
        '''
        return 2 * len(self.recognizers) - 1
    
    def count_connectors(self) -> int:
        '''
        Returns the number of connectors of the organism
        '''
        return len(self.connectors)

    def count_recognizers(self) -> int:
        '''
        Returns the number of recognizers of the organism
        '''
        return len(self.recognizers)
    
    def sum_pssm_lengths(self) -> int:
        '''
        Returns the sum of the lengths of all the PSSMs of the organism.
        '''
        sum_lengths = 0
        for pssm in self.recognizers:
            sum_lengths += pssm.length
        return sum_lengths
    
    def get_gap_score(self, connector_idx, d, s_dna_len):
        """Calls the appropriate connector, with the given distance and length of the DNA sequence to
		   obtain the energy of the connector.
		"""
        if d < s_dna_len:
            # !!! New input param for connector.get_score: the list of the lengths of the recognizers
            recog_lengths = [recog.length for recog in self.recognizers]
            gap_score = self.connectors[connector_idx].get_score(d, s_dna_len, recog_lengths)
            return gap_score
        else:
            return -1 * np.inf
    
    def get_score_from_pssm(self, row_idx_from_placement_matrix, nucleotide):
        """Calls the appropriate PSSM (and column) to obtain the score, given a nucleotide
		"""
        pssm_index = self.row_to_pssm[row_idx_from_placement_matrix][0]
        pssm_column = self.row_to_pssm[row_idx_from_placement_matrix][1]
        pssm_object = self.recognizers[pssm_index]
        score = pssm_object.pssm[pssm_column][nucleotide]
        return score
    
    def get_diag_score(self, pointers_mat, row_idx, col_idx, dna_sequence):
        """Evaluates and returns a substitution score (diagonal move), using
		   get_score_from_pssm and taking into account several special cases.
		   row_idx and col_idx identify the "destination" cell [the cell being
		   evaluated]
		"""
    
        diag_score = 0
        
        # Check if it is a gap of zero bp
		# This means two PSSMs back to back, which then must incorporate a zero
		# gap score
        if self.is_a_0_bp_gap(pointers_mat, row_idx, col_idx):
            # Call connector, for a zero bp gap evaluation, add it to the
			# diagonal score [connector needs to the length of the DNA seq]
            pssm_idx = self.row_to_pssm[row_idx][0]
            connector = self.connectors[pssm_idx - 1]
            
            # !!! New input param for connector.get_score: the list of the lengths of the recognizers
            recog_lengths = [recog.length for recog in self.recognizers]
            
            zero_gap_score = connector.get_score(0, len(dna_sequence), recog_lengths)
            diag_score += zero_gap_score
        
		# get nucleotide and compute PSSM score for it
        nucleotide = dna_sequence[col_idx - 1] 
        pssm_score = self.get_score_from_pssm(row_idx, nucleotide)
        diag_score += pssm_score
        
        return diag_score
    
    def is_first(self, row_idx_from_placement_matrix):
        """Returns True if we are on the first element of a PSSM recognizer
		"""
        pssm_col = self.row_to_pssm[row_idx_from_placement_matrix][1]
        if pssm_col == 0:
            return True
        else:
            return False
    
    def is_last(self, row_idx_from_placement_matrix): 
        """Returns true if we are on the last element of a PSSM recognizer
		"""
		# if next one is a first, then we are at a last ;-)
        if self.is_first(row_idx_from_placement_matrix + 1):
            return True
        else:
            return False
    
    def is_a_0_bp_gap(self, pointers_mat, row_idx, col_idx):
        """Tells whether the cell defines a contiguous diagonal
		   run between two PSSMs
		"""

		# if the row does not correspond to the first column of a PSSM
        if self.is_first(row_idx)==False:
            return False
        
		# if the row IS a first column of a PSSM
		# get row index of diagonal-up-left cell
		# this should be the row of the cell pointing to the
		# end of the previous PSSM
        pointer_row_idx = pointers_mat[0][row_idx-1, col_idx-1]
        
		# the equality below, will only be true, if there was
		# a diagonal move. this, combined with the fact that we know
		# that this was a PSSM last row, identifies that the diagonal
		# up-left element comes from a PSSM diagonal score
		# (so you really have a back-to-back PSSM situation)
        if pointer_row_idx == row_idx-2:
            return True
        else:
            return False
    
    def get_random_connector(self) -> int:
        """Returns the index of a random connector of the organism

        Returns:
            Integer between 0 and N-1 (both included), where N is the number of
            connectors the organism has.
        """
        
        num_connectors =  self.count_connectors()
        return random.randint(0, num_connectors - 1)


    def get_random_recognizer(self) -> int:
        """Returns the index of a random recognizer of the organism

        Returns:
            Integer between 0 and N-1 (both included), where N is the number of
            recognizers the organism has.
        """
        
        num_recognizers =  self.count_recognizers()
        
        return random.randint(0, num_recognizers - 1)
    
    def set_connectors(self, connectors_list):
        """Set the connectors of the organism to be those provided in the input
        list.
        """
        self.connectors = connectors_list
        
    def set_recognizers(self, recognizers_list):
        """Set the recognizers of the organism to be those provided in the
        input list.
        """
        self.recognizers = recognizers_list
        # Set/update the row_to_pssm attribute, used for the placement
        self.set_row_to_pssm()
    
    def append_connector(self, connector_obj):
        """ Adds a copy of the given connector to the organism, by appending
        it to the  organism.connectors  list.
        """
        connector = copy.deepcopy(connector_obj)
        self.connectors.append(connector)

    def append_recognizer(self, recognizer_obj):
        """ Adds a copy of the given recognizer to the organism, by appending
        it to the  organism.recognizers  list.
        """
        recognizer = copy.deepcopy(recognizer_obj)
        self.recognizers.append(recognizer)
        # Set/update the row_to_pssm attribute, used for the placement
        self.set_row_to_pssm()
        
    def print(self) -> None:
        """Prints the whole tree data structure
        """
        
        print("***** Organism {} *****".format(self._id))
        for i in range(len(self.recognizers) - 1):
            self.recognizers[i].print()
            self.connectors[i].print()
        self.recognizers[-1].print()

    def export(self, filename: str) -> None:
        """Exports the whole tree data structure

        Args:
            filename: Name of the file to export the organism
        """
        organism_file = open(filename, "a+")
        organism_file.write("***** Organism {} *****".format(self._id))
        
        for i in range(len(self.recognizers) - 1):
            self.recognizers[i].export(organism_file)
            self.connectors[i].export(organism_file)
        self.recognizers[-1].export(organism_file)

        organism_file.write("\n\n")
        organism_file.close()

    def export_results(self, a_dna: list, filename: str) -> None:
        """Exports the binding profile of the organism against each of the 
           DNA sequences provided as a list

        Args:
            filename: Name of the file to export sequences
            a_dna: list fo sequences to export

        """
        
        ofile = open(filename, "a+")
        # for each DNA sequence
        for s_dna in a_dna:
            placement = self.get_placement(s_dna)
            placement.print_placement(outfile = ofile)
        ofile.close()

    def print_result(self, s_dna: str) -> None:
        """Prints the binding profile of the organism against the 
           provided DNA sequence 
           
        Args:
            s_dna: DNA sequence to export

        Returns:
            DNA sequence and binding sites of the organisms recognizer
        """

        s_dna = s_dna.lower()
        placement = self.get_placement(s_dna)
        placement.print_placement(stdout = True)

    def adjust_connector_scores(self, c_idx, c_scores, sequence_length):
        #get to the start in the array of connector info
        offset = sum([(self.connectors[i].expected_seq_length * 2 + 2) for i in range(0, c_idx)])

        #get to first precomputed position
        offset += 2

        #get to the length of the sequence
        offset += sequence_length

        con = self.connectors[i]
        g_curr = g(sequence_length, con._mu, con._sigma)
        h_curr = 1

        h_list = [1 + con.pseudo_count]
        tot_sum =  h_list[0]

        for i in range(sequence_length - 1, -1, -1):
            g_next = g(i, mu, sigma)
            h_next = h_curr * math.exp(g_next - g_curr)

            h_list.append(h_next + con.pseudo_count)
            tot_sum += h_next + con.pseudo_count

            g_curr = g_next
            h_curr = h_next

        h_list_sorted = h_list[::-1]
        p_list = [h/tot_sum for h in h_list_sorted]

        prob_sum = 0.0
        for i in range(len(p_list)):
            prob_sum += p_list[i]
            c_scores[offset + i] = np.log(p_list[i])
            c_scores[offset + con.expected_seq_length + i] = np.log(prob_sum)

    def get_placement(self, sequence: str) -> PlacementObject:

        """
        Calculates a placement for a given organism for a given sequence using
        the _calculatePlacement module.

        Args:
            sequence: DNA sequence to place organism on

        Returns:
            PlacementObject containing information of optimal placement
        """
        if self.minimum_length > len(sequence):
            placement = Placement(self._id, sequence)
            placement.set_energy(-1E100)
            return placement

        number_PSSM = len(self.recognizers)
        max_length = -1
        if number_PSSM > 1:
            max_length = self.connectors[0].expected_seq_length - 1
            for i in range(len(self.connectors)):
                if len(sequence) < self.connectors[i].adjust_score_threshold:
                   self.adjust_connector_scores(i, len(sequence))
                     
        # Get an array of lengths of each recognizer in the organism

        # instantiation of numpy arrays that will hold placement info
        # gathered by the _calculatePlacement module
            
        
        gaps = np.empty(number_PSSM, dtype = np.dtype('i'))
        gap_scores = np.empty(number_PSSM - 1, dtype = np.dtype('d'))
        PSSM_scores = np.empty(number_PSSM + 1, dtype = np.dtype('d'))

        c_scores = np.copy(self.connectors_scores_flat)
        if number_PSSM > 1:
            max_length = self.connectors[0].expected_seq_length - 1
            for i in range(len(self.connectors)):
                if len(sequence) < self.connectors[i].adjust_score_threshold:
                   self.adjust_connector_scores(i, c_scores)
            
        _multiplacement.calculate(bytes(sequence, "ASCII"), bytes(self.recognizer_types, "ASCII"), self.recognizers_flat, self.recognizer_lengths,  c_scores, PSSM_scores, gap_scores, gaps, max_length, self.recognizer_models, self.recognizer_bin_edges, self.recognizer_bin_nums)

        # parse data from the _calculatePlacement module and put it
        # into a PlacementObject to be returned
        placement = PlacementObject(self._id, sequence)
        placement.set_energy(float(PSSM_scores[-1]))
        placement.set_recognizer_scores([float(score) for score in PSSM_scores[:-1]])
        placement.set_connectors_scores([float(score) for score in gap_scores])

        current_position = gaps[0]
        placement.set_recognizer_types(self.recognizer_types)
        for i in range(number_PSSM):

            stop = current_position + self.recognizer_lengths[i]
            placement.append_recognizer_position([int(current_position), int(stop)])
            current_position += self.recognizer_lengths[i]

            if i < number_PSSM - 1:
                stop = current_position + gaps[i + 1]
                placement.append_connector_position([int(current_position), int(stop)])
                current_position += gaps[i + 1]
        
        #placement.print_placement(stdout=True)
        #print("python overhead took {} seconds".format((t2 + t1) *10E-9))
        """
        for i in ['m', 't', 'h', 'r']:
            if i in self.recognizer_types:
                placement.print_placement(stdout=True)
        """
        return placement

    def set_minimum_length(self) -> None:
        self.minimum_length = sum([i.length for i in self.recognizers])

    def set_pf_auc(self) -> None:
        s = time.monotonic_ns()
        #if we have no connectors, empty connectors info
        n_con = len(self.connectors)
        if n_con == 0:
            self.connectors_scores_flat = np.zeros(0, dtype=np.dtype('d'))
            return

        #allocate all of the space needed for the whole flat array
        exp_len  = self.connectors[0].expected_seq_length
        self.connectors_scores_flat = np.zeros(n_con * (2 * exp_len + 2), dtype=np.dtype('d'))
        offset = 0
        pf = 0
        auc = 0

        #this value keeps track of the the length n that a sequence must be
        #such that all values from index 0 to index n in the precomputed scores
        #are all rounded to con.pseudo_count
        adjust_score_threshold = 0
        """
        for each connector we need to:
            add mu and sigma
            compute the cdf(-0.5) for normalization
        """

        for i in range(n_con):
            con = self.connectors[i]
            prev_cdf = norm_cdf(-0.5, con._mu, con._sigma)
            cdf_0 = prev_cdf
            self.connectors_scores_flat[offset] = np.double(con._mu)
            self.connectors_scores_flat[offset + 1] = np.double(con._sigma)
            offset += 2

            """
            range is 1 to e_len + 1 because first pdf is cdf(0.5) - cdf(-0.5)
            and the last pf should be cdf(e_len + 0.5) - cdf(e_len - 0.5)
            since we're starting at 1, for each index in the auc array we need
            to actually use the previous iteration of the auc since the first index
            should be cdf(-0.5) - cdf(-0.5)
            """

            for j in range(1, exp_len + 1):
                cdf = norm_cdf(j - 0.5, con._mu, con._sigma)
                if cdf == prev_cdf and cdf <= con.pseudo_count:
                    con.adjust_score_threshold = j
                auc = np.log2(prev_cdf - cdf_0 + (con.pseudo_count * j))

                pf = np.log2(cdf - prev_cdf + con.pseudo_count)

                self.connectors_scores_flat[offset] = np.double(pf)
                self.connectors_scores_flat[offset + exp_len] = np.double(auc)
                offset += 1
                prev_cdf = cdf

            #our offset to get to the start of the next connector will be e_len since
            #we have alread iterated e_len times and each connector takes (2 * e_len + 2)
            #indices
            offset += exp_len
        #print("setting stuff took:", (time.monotonic_ns() - s) * 1E-9)
        
    def flatten(self) -> None:
        """
        Flattens an organisms recognizers into a 1D array
        computes gap scores and forms a 1D array out of them.
        This function is to be called whenever an organism is made
        or changed.
        
        Args:
            None
        
        Returns:
            None
        """
        import time
        # instantiation of lists to hold flattened info
        s = time.monotonic_ns()
        flat_recognizers = []
        flat_connector_scores = []
        flat_rec_models = []
        rec_bin_nums = []
        rec_bin_edges = []
        recognizer_lengths = []
        self.recognizer_types = ""
        
        for recognizer in self.recognizers:
            if recognizer.get_type() == 'p':
                for column in recognizer.pssm:
                    for base in ['a','g','c','t']:
                        flat_recognizers.append(column[base])
            else:
                for prob in recognizer.null_model:
                    flat_rec_models.append(prob)
                for prob in recognizer.alt_model:
                    flat_rec_models.append(prob)
                for edge in recognizer.bins:
                    rec_bin_edges.append(edge)
                rec_bin_nums.append(len(recognizer.bins))
            recognizer_lengths.append(recognizer.length)
            self.recognizer_types += recognizer.get_type()
        # organism holds a numpy array of the flattened lists
        self.recognizers_flat = np.array(flat_recognizers, dtype = np.dtype('d')) 
        self.recognizer_lengths = np.array(recognizer_lengths, dtype = np.dtype('i'))
        self.recognizer_models = np.array(flat_rec_models, dtype = np.dtype('d'))
        self.recognizer_bin_nums = np.array(rec_bin_nums, dtype = np.dtype('i'))
        self.recognizer_bin_edges = np.array(rec_bin_edges, dtype = np.dtype('d'))

        self.set_pf_auc()
        return

        """
        g_s = time.monotonic_ns()
        if self.is_precomputed == True:
            for connector in self.connectors:
                flat_connector_scores += [connector._mu, connector._sigma]
                flat_connector_scores += (connector.stored_pdfs + connector.stored_cdfs)
#               for i in range(connector.expected_seq_length):

            # organism holds a numpy array of the flattened lists
            self.connectors_scores_flat = np.array(flat_connector_scores, dtype = np.dtype('d'))

        else:
            con_mu_sigma = []
            for connector in self.connectors:
                con_mu_sigma.append(connector._mu)
                con_mu_sigma.append(connector._sigma)

            self.connectors_scores_flat = np.array(con_mu_sigma, dtype=np.dtype('d'))
        e = time.monotonic_ns()
        print("full gap precomputation and population took:", (e-g_s) * 10E-9)
        #print("flattening organism took {} seconds".format((e - s) *10E-9))
        """