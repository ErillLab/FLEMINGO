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
    """
    Organism object
    The Organism object essentially contains two vectors:
        - A vector of connector objects
        - A vector of recognizer objects
    The order of the elements in these vectors determine, implicitly, the
    connections between the elements.
    """
    
    def __init__(self, _id: str, conf: dict) -> None:
        """Organism constructor

        Args:
            _id: Organism identifier (assigned by factory)
            conf: Organism-specific configuration values from JSON file
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
        self.sum_recognizer_lengths = 0
        
        self.is_precomputed = conf["PRECOMPUTE"]
	
        #probability of deleting/inserting a recognizer
        self.probability_delete_recognizer = conf["PROBABILITY_DELETE_RECOGNIZER"]
        self.probability_insert_recognizer = conf["PROBABILITY_INSERT_RECOGNIZER"]
	
        # probability of mutating a node
        self.probability_node_mutation = conf["PROBABILITY_NODE_MUTATION"]
        
        # probability of applying a nudge shift to a recognizer
        self.probability_nudge_recognizer = conf["PROBABILITY_NUDGE_RECOGNIZER"]
        
        # min and maximum number of recognizers allowed
        self.min_num_of_recognizers = conf["MIN_NUM_OF_RECOGNIZERS"]
        self.max_num_of_recognizers = conf["MAX_NUM_OF_RECOGNIZERS"]
        
        # !!! maximum length of PSSMs allowed
        # self.max_pssm_length = max_pssm_length
        
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
                                      'connectors': None,
                                      'connectors_adjustments': None}
        
        self.fitness = None
    
    # def set_assembly_instructions_old_version(self, aligned_repres, connectors_table, p1_id, p2_id):
    #     '''
    #     Sets the self.assembly_instructions attribute.
    #     '''
    #     self.assembly_instructions['p1'] = p1_id
    #     self.assembly_instructions['p2'] = p2_id
    #     self.set_recogs_assembly_instructions(aligned_repres)
    #     self.set_connectors_assembly_instructions(aligned_repres, connectors_table)
    #     self.set_connectors_adjustments(aligned_repres)
    
    def set_assembly_instructions(self, aligned_repres, organism_tag):
        '''
        !!! Docstring here ...
        '''
        # Parent IDs
        self.assembly_instructions['p1'] = aligned_repres.parents_ids[0]
        self.assembly_instructions['p2'] = aligned_repres.parents_ids[1]
        
        # Child representation
        if organism_tag == 'org1':
            child_repres = aligned_repres.organism1
        elif organism_tag == 'org2':
            child_repres = aligned_repres.organism2
        else:
            raise ValueError("Unknown organism_tag. It should be 'org1' or 'org2'.")
        
        self.set_recogs_assembly_instructions(child_repres)
        self.set_connectors_assembly_instructions(child_repres, aligned_repres.connectors_table)
        self.set_connectors_adjustments(aligned_repres)
    
    def set_recogs_assembly_instructions(self, child_repres):
        '''
        Compiles the list of recognizers in the self.assembly_instructions
        attribute.
        '''
        child_recog_names = [item for item in child_repres if item != '-']
        recogs_list = []
        for item in child_recog_names:
            if item not in recogs_list:
                recogs_list.append(item)
        self.assembly_instructions['recognizers'] = recogs_list
    
    def set_connectors_assembly_instructions(self, child_repres, connectors_table):
        '''
        Compiles the list of connectors in the self.assembly_instructions
        attribute.
        '''
        required_connectors = self.annotate_required_connectors(child_repres)
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
    
    def set_connectors_adjustments(self, aligned_repres):
        '''
        !!! Docstring here ...
        '''
        self.assembly_instructions['connectors_adjustments'] = {
            'p1_connectors': aligned_repres.org1_connectors_adjustments,
            'p2_connectors': aligned_repres.org2_connectors_adjustments}
    
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
        '''
        The `row_to_pssm` attribute maps each row of the alignment
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
        '''
        
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

    def set_id(self, _id: str) -> None:
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
        if random.random() < self.probability_delete_recognizer:
            # Delete a recognizer (and a connector)
            self.delete_recognizer()
        
        # Insertion
        if random.random() < self.probability_insert_recognizer:
            # Insert a recognizer (and a connector)
            self.insert_recognizer(org_factory)
        
        # Nudge
        if random.random() < self.probability_nudge_recognizer:
            # Nudge a recognizer (two connectors are modified)
            self.nudge_recognizer()
        
        # Node-mutations
        self.mutate_nodes(org_factory)
        
        # After applying mutations, order/columns/number of pssm's may have
        # changed, so we call the set_row_to_pssm to set their mapping on the
        # alignment matrix anew
        self.set_row_to_pssm()
        self.flatten()
    
    def nudge_recognizer(self, shift='random'):
        '''
        When the organism is placed on a sequence, the location of a recognizer
        in the optimal placement depends on its binding energy to the different
        DNA k-mers available, but also on the pulling/pushing forces from the
        adjacent connectors. This function can be used to shift the expected
        location of a recognizer, relative to the previous one and the next one.
        Unless specified otherwise, the shift is a random number between 0 and
        1, i.e., it's a small nudge.
        
        Args:
            shift : optional
                The default is "random". When it's "random" (not specified) the
                shift amount is a random number between 0 and 1 (`shift` is
                redefined). This is equivalent to a small nudge. Otherwise, a
                numeric value can be specified for the `shift`.
        
        Example:
            
            Consider the following placement
            
            
            CTGT........ACAG...............GTTC
            pppp--------pppp---------------pppp
            _____mu=8.2__________mu=14.8_______
            
            
            The average distance between the first and the third recognizer is
            optimal, and it is 8.2 + 4 + 14.8 = 27. However, it would be more
            fit if it was 8 + 4 + 15 = 27.
            That's what we would get if we apply this function to the middle
            recognizer with a `shift` of -0.2.
        '''
        # Choose a recognizer
        idx = random.randint(0, self.count_recognizers() - 1)
        
        # A small shift (between 0 and 1) is applied, unless specified otherwise
        if shift == 'random':
            shift = random.random()
        
        # Decide if the recognizer is puhsed to the left or to the right
        # Negative shift means shift towards the left
        if random.random() < 0.5:
            shift = -1 * shift
        
        # If it's not the first recognizer, the connector to the left exists,
        # and it needs to be mutated
        if idx > 0:
            # Apply the shift to the connector to the left
            self.connectors[idx-1]._mu += shift
            self.connectors[idx-1].set_precomputed_pdfs_cdfs()
        
        # If it's not the last recognizer, the connector to the right exists,
        # and it needs to be mutated
        if idx < self.count_recognizers() - 1:
            # Apply the inverse shift to the connector to the right
            # (if the connector to the left shrinks, the one to the right
            # expands, and vice versa)
            self.connectors[idx]._mu -= shift
            self.connectors[idx].set_precomputed_pdfs_cdfs()
    
    def delete_recognizer(self, recognizer_to_remove=None):
        '''
        A randomly selected recognizer is deleted, unless the index of a
        specific target recognizer to be deleted is provided (the
        `recognizer_to_remove` argument).
        
        Connectors are adjusted as needed:
        ----------------------------------
            
            - If the deleted recognizer was a terminal node, its adjacent
              connector is also deleted.
              
            - If the deleted recognizer was an internal node, its two adjacent
              connectors are replaced by a single connector with mu and sigma
              such that it models the spacer equivalently to the sum of the two
              pre-existing connectors together with the deleted recognizer.
        
        Args:
            recognizer_to_remove : int, optional
                Index of the recognizer to be deleted. The default is None.
                When None (not specified) a random recognizer is deleted.
        '''
        n_recognizers = self.count_recognizers()
        
        # if this is a single-recognizer organism, no deletion occurs
        if n_recognizers < 2:
            return None
        
        # Unless specified, choose randomly the recognizer to be deleted
        if recognizer_to_remove is None:
            recognizer_to_remove = random.randint(0, n_recognizers - 1)
        
        connector_to_keep = None
        
        # Establish what connector should be removed
        if recognizer_to_remove == 0:
            # If the recognizer to delete is the first, the connector to
            # delete is the one to the right
            # ( same index: connector_to_remove = recognizer_to_remove = 0 )
            connector_to_remove = recognizer_to_remove
        elif recognizer_to_remove == n_recognizers - 1:
            # If the recognizer to delete is the last, the connector to
            # delete is the one to the left
            # ( different index: connector_to_remove = recognizer_to_remove - 1 )
            connector_to_remove = recognizer_to_remove - 1
        else:
            # if the recognizer to delete is in not a terminal recognizer of
            # the chain, the parent connector to be deleted (left/right) is
            # chosen randomly
            if random.random() < 0.5:
                # Index of the parent connector that will be deleted
                connector_to_remove = recognizer_to_remove
                # Index of the parent connector that will be adjusted
                connector_to_keep = recognizer_to_remove - 1
            else:
                # Index of the parent connector that will be deleted
                connector_to_remove = recognizer_to_remove - 1
                # Index of the parent connector that will be adjusted
                connector_to_keep = recognizer_to_remove
        
        # "intelligent" method: a new connector is created by merging the
        # connectors on the sides of the deleted recognizer into a single
        # larger connector. In practice, the connector that was not removed
        # (connector_to_keep) is expanded to cover the span that was covered by
        # the two connectors (connector_to_remove and connector_to_keep) + the
        # length of the deleted recognizer
        
        if connector_to_keep:
            # Adjust parameters of the neighbour connector
            ''' The parent connector that is not deleted is modified,
            so that it can span the gap left by the deletion, without
            heavily affecting the placement of the nodes to the sides of
            the deletion point.'''
            
            # Adjust MU
            '''Thus, we need to add to the mu of the remaining connector
            also the mu of the deleted connector and the legth of the
            deleted recognizer.'''
            adj_mu = (self.connectors[connector_to_keep]._mu +
                      self.connectors[connector_to_remove]._mu +
                      self.recognizers[recognizer_to_remove].length)

            
            # Adjust SIGMA
            ''' Variances are additive (under assumption of independence
            between the two random variables). Therefore to adjust the
            variance we need to add the variance of the deleted connector
            (the length of the deleted recognizer has no variance). Therefore,
            its standard deviation (sigma) is the square root of the sum of the
            squares of the sigmas. '''
            adj_sigma = (self.connectors[connector_to_keep]._sigma ** 2 +
                        self.connectors[connector_to_remove]._sigma ** 2)**(1/2)
            
            # set new mu and new sigma
            self.connectors[connector_to_keep].set_mu(adj_mu)
            self.connectors[connector_to_keep].set_sigma(adj_sigma)
            self.connectors[connector_to_keep].set_precomputed_pdfs_cdfs()
        
        # Recreate arrays of recognizers and connectors skipping the recgonizer
        # and the connector selected for deletion				
        new_recognizers = (self.recognizers[:recognizer_to_remove] +
                           self.recognizers[recognizer_to_remove + 1:])
        new_connectors =  (self.connectors[:connector_to_remove] +
                           self.connectors[connector_to_remove + 1:])
        
        # assign new arrays
        self.recognizers = new_recognizers
        self.connectors = new_connectors
    
    def insert_recognizer(self, org_factory):
        '''
        A new recognizer is inserted in a random position of the chain.
        
        Connectors are adjusted as needed:
        ----------------------------------
            
            - If the inserted recognizer is a terminal node, a new random
              connector is used to link it to the chain.
              
            - If the recognizer is inserted internally in the chain, the
              connector modeling the spacer where the insertion happened is
              replaced by two connectors such that, together with the newly
              inserted recognizer, they span the same spacer (and with the same
              variablity) as the one specified by the pre-existing connector.
        
        Args:
            org_factory : OrganismFactory
        '''
        
        n_recognizers = self.count_recognizers()
        
        # instantiate the new recognizer and the new connector
        new_connector = org_factory.create_connector()
        new_recognizer = org_factory.create_recognizer()
        
        # Special case for empty organism.
        # This can happen after a recombination event where all the parental
        # nodes were assigned to the other child, leaving this child empty
        if n_recognizers == 0:
            # In this case, the connector must not be added
            self.recognizers.append(new_recognizer)
            return None
        
        # Choose randomly the recognizer next to which the insertion
        # is going to occur
        recognizer_idx = random.randint(0, n_recognizers - 1)
        
        
        # "intelligent" method: the pre-existing connector and the new one are
        # "shrinked" so that together they model a spacer "equivalent" to
        # the one modeled by the connector previously occupying the space
        # where the new recognizer has been inserted
        
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
            self.connectors[connector_to_compress].set_precomputed_pdfs_cdfs()
            # Update paramters of inserted connector
            new_connector.set_mu(inserted[0])
            new_connector.set_sigma(inserted[1])
            new_connector.set_precomputed_pdfs_cdfs()
        
        # recreate arrays of connectors and recognizers, adding
        # the newly minted recognizer+connector and comprising
        # also the "compressed" pre-existing connector
        new_recognizers = (self.recognizers[:recognizer_idx] +
                           [new_recognizer] +
                           self.recognizers[recognizer_idx:])
        new_connectors = (self.connectors[:connector_idx] +
                           [new_connector] +
                           self.connectors[connector_idx:])
        
        # assign new connectors and recognizers arrays
        self.recognizers = new_recognizers
        self.connectors = new_connectors
    
    def mutate_nodes(self, org_factory):
        '''
        Triggers node-specific mutations on the nodes of the organism.
        
        If `probability_node_mutation` is specified, only one node will
        become prone to mutations, with a probability equal to
        `probability_node_mutation`. Otherwise, all the nodes are always
        prone to mutations.
        
        When a node is prone to mutations it doesn't mean it will necessarily
        be mutated. Each type of node mutation can occur with its own specific
        proability (specified in the config).
        '''
        
        # If PROBABILITY_NODE_MUTATION is set to real value, a single
        # node is selected. If set to null, all nodes are prone to mutations
        if self.probability_node_mutation:
            # Mutate a random node
            if random.random() < self.probability_node_mutation:
                
                n_nodes = self.count_nodes()
                random_node_idx = random.randint(0, n_nodes - 1)
                if random_node_idx < self.count_recognizers():
                    # mutate a recognizer
                    moved_pssm_bounds = self.recognizers[random_node_idx].mutate(org_factory)
                    # Adjust adjacent gaps if needed
                    if moved_pssm_bounds:
                        self.adjust_gaps_after_recog_bounds_displacement(
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
                    self.adjust_gaps_after_recog_bounds_displacement(i, moved_pssm_bounds)
            
            # Connectors
            for connector in self.connectors:
                connector.mutate(org_factory)

    def adjust_gaps_after_recog_bounds_displacement(self, recog_idx, displacement_code):
        '''
        For well-adapted organisms (organisms that have a consistent placement
        strategy on the positive set) some recognizer mutations will produce a
        change in the DNA positions occupied by the recognizers.
        Example: If a PSSM gets a new column added to the left, and if there's
        a gap to the left of that PSSM, the connector modeling that gap should
        decrease its mu by 1, so that it preserves its optimality. Otherwise,
        the benefit of a good extra PSSM column could be overcomed by the
        downsides of having the average gap size off by 1 unit in the adjacent
        connector.
        
        The recognizer mutations to be taken into account for this problem are:
            - increase/decrease PSSM size
            - increase/decrease SHAPE recognizer size
            - right/left shift of the PSSM
        Those are the mutation that can produce some displacement of the
        recognizer's bounds (the left and the right bounds). When they're
        called, a code is generated, that keeps track of how the bounds have
        changed. This code is the second input argument of this function. The
        purpose of this function is to adjust connectors according to those
        codes, so that recognizers are free to change size or shift without
        damaging the optimality of adjacent connectors.

        Parameters
        ----------
        recog_idx : int
            Index of the recognizer that has gone through bounds displacement.
        displacement_code : list
            It's a list of two integers. The first integer is the shift
            occurred the left bound of the recognizer. The second integer is
            the shift occurred to its right bound.
            A positive shift means that the bound moved to the right, while a
            negative shift means a shift to the left.
            The `displacement_code` is generated when a recognizer is mutated
            by its `mutate()` method.
        '''
        
        left_bound_shift, right_bound_shift = displacement_code
        
        # If the LEFT bound has changed
        if left_bound_shift != 0:
            # Adjust connector to the left
            conn_index = recog_idx - 1
            
            # There's no connector to the left of the first recognizer
            # (no connector with index -1), so check that the index is not -1.
            if conn_index >= 0:
                # Adjust mu
                new_mu = self.connectors[conn_index]._mu + left_bound_shift
                self.connectors[conn_index].set_mu(max(0, new_mu))
                
                # Update connector's PDF and CDF values
                self.connectors[conn_index].set_precomputed_pdfs_cdfs()
        
        # If the RIGHT bound has changed
        if right_bound_shift != 0:
            # Adjust connector to the right
            conn_index = recog_idx
            
            # There's no connector to the right of the last recognizer
            # so check that the index is not equal to the number of connectors
            if conn_index < self.count_connectors():
                # Adjust mu
                new_mu = self.connectors[conn_index]._mu - right_bound_shift
                self.connectors[conn_index].set_mu(max(0, new_mu))
                
                # Update connector's PDF and CDF values
                self.connectors[conn_index].set_precomputed_pdfs_cdfs()
        self.flatten()
    
    def check_size(self, org_factory):
        '''
        Checks that the organism complies with the constraints specified in the
        config file. If the organism doesn't comply, instead of killing it (which
        means that the parent it's paired with doesn't have a competitor), the
        organism is modified as necessary to make it acceptable.
        
        Constraints:
            - the number of nodes in the organism
            - the sum of the lengths of the recognizers
        
        '''
        
        edited = False
        
        # If the number of nodes exceeded the maximum allowed, prune the organism
        if self.max_num_of_recognizers != None:
            while self.count_recognizers() > self.max_num_of_recognizers:
                self.delete_recognizer()
                edited = True
        
        # If the number of nodes is not sufficient, insert recognizers
        while self.count_recognizers() < self.min_num_of_recognizers:
            self.insert_recognizer(org_factory)
            edited = True
        
        # If a connector has mu larger than the max possible gap, reduce mu
        # (such large mu values can be generated by the deletion operator)
        for conn in self.connectors:
            if conn._mu > org_factory.max_seq_length - 2:
                conn._mu = org_factory.max_seq_length - 2
                conn.set_precomputed_pdfs_cdfs()
                edited = True
        
        # If the recognizers cover more bp than available on the shortest sequence,
        # shrink a random recognizer
        while sum([r.length for r in self.recognizers]) > org_factory.min_seq_length:
            # Choose a random recognizer and decrease its length
            idx = random.randint(0, self.count_recognizers()-1)
            displacement_code = self.recognizers[idx]._decrease_length([0, 0])
            if displacement_code:
                self.adjust_gaps_after_recog_bounds_displacement(idx, displacement_code)
            # If it was not possible to decrease the length (it was already
            # at the minimum length allowed for a recognizer) delete it
            else:
                self.delete_recognizer(recognizer_to_remove=idx)
            edited = True
        
        # Update organism
        if edited:
            self.set_row_to_pssm()
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
        also make sense under the assumption that discrimination happens in a
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
            gamma*100 is the percentage of discarded elements from the left tail
            of the distribution
        Trim-left-M-SD:
            TRIMMED MEAN AND TRIMMED STANDARD ERROR.
            gamma*100 is the percentage of discarded elements from the left tail
            of the distribution
        
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
        
        if n_to_trim == 0 or method == "Welch":
            # Mean
            mean = np.mean(energy_scores)
            # Standard deviation
            stdev = np.std(energy_scores, ddof=1)
            stdev = max(1, stdev)  # Lower bound to sigma (see docstring)
            sterr = stdev / np.sqrt(len(energy_scores))
        elif method == "Yuen":
            # Trimmed mean
            # Make sure that `n_to_trim` doesn't exceed the maximum. We should always
            # trim less than half of the elements from one side, to avoid being left
            # with an empty list
            max_n_to_trim = int((len(energy_scores)-1)/2)
            n_to_trim = min(n_to_trim, max_n_to_trim)
            energy_scores.sort()
            energy_scores_trimmed = energy_scores[n_to_trim:-n_to_trim]
            mean = np.mean(energy_scores_trimmed)
            # Standard deviation
            stdev = np.std(energy_scores, ddof=1)
            stdev = max(1, stdev)  # Lower bound to sigma (see docstring)
            sterr = stdev / np.sqrt(len(energy_scores))
        elif method == "Trim-left-M":
            # Left-trimmed mean
            energy_scores.sort()
            energy_scores_trimmed = energy_scores[n_to_trim:]
            mean = np.mean(energy_scores_trimmed)
            # Standard deviation
            stdev = np.std(energy_scores, ddof=1)
            stdev = max(1, stdev)  # Lower bound to sigma (see docstring)
            sterr = stdev / np.sqrt(len(energy_scores))
        elif method == "Trim-left-M-SD":
            # Left-trimmed mean
            energy_scores.sort()
            energy_scores_trimmed = energy_scores[n_to_trim:]
            mean = np.mean(energy_scores_trimmed)
            # Left-trimmed standard deviation
            stdev = np.std(energy_scores_trimmed, ddof=1)
            stdev = max(1, stdev)  # Lower bound to sigma (see docstring)
            sterr = stdev / np.sqrt(len(energy_scores_trimmed))
        else:
            raise ValueError('Unknown method "' + method + '". Please check ' +
                             'the FITNESS_FUNCTION paramter in the config file.')
        return mean, sterr
    
    def set_fitness(self, pos_set, neg_set, method, gamma=0.2):
        '''
        Sets the fitness of the organism.

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
        self.fitness = (M_p - M_n) / (SE_p**2 + SE_n**2)**(1/2)
    
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
        placement = self.get_placement(s_dna)
        placement.print_placement(stdout = True)

    def adjust_connector_scores(self, c_idx, c_scores, sequence_length):
        #####DEPRECATED####
        """
        Rescales the probabilties for a connector when all possible gap lengths are
        below a certain threshold (this could cause the placement function to 
        improperly treat all gap lengths as having the same probability in the 
        alternative model)

        Args:
            c_idx: the index of the connector in self.connectors that we wish to rescael
            c_scores: the 1D array of mus, sigmas, pfs, and aucs for all connectors
            sequence_length: the length of the sequence that we wish to place the organism on
        """
        # con is just used as a shorthand to refer to the desired connector
        con = self.connectors[c_idx]
        #offset = sum([(self.connectors[c_idx].max_seq_length * 2 + 2) for i in range(0, c_idx)]) + 2
        #self.connectors[c_idx].adjust_scores(c_scores[offset:], sequence_length, self.sum_recognizer_lengths)
        #return
        # treat sigma to be 1E-100 if it is lower, otherwise overflow potential
        # is greater
        sigma = min(con._sigma, 1E-100)
        g_curr = g(sequence_length - self.sum_recognizer_lengths, con._mu, sigma)
        h_curr = 1

        h_list = [1 + con.pseudo_count]
        tot_sum =  h_list[0]

        for i in range(sequence_length - 1 - self.sum_recognizer_lengths, -1, -1):
            g_next = g(i, con._mu, sigma)
            h_next = h_curr * math.exp(g_next - g_curr)

            h_list.append(h_next + con.pseudo_count)
            tot_sum += h_next + con.pseudo_count

            g_curr = g_next
            h_curr = h_next

        h_list_sorted = h_list[::-1]
        p_list = [h/tot_sum for h in h_list_sorted]

        # get to the start in the array of connector info
        # each connector occupies indices for its mu, sigma, pfs, and aucs
        # this means each connector uses 2 + (indices for pfs) + (indices for aucs)
        # which is 2 + 2 * max_seq_length
        # this computation is done for each connector before the connector we need to re-scale in order
        # to get the starting index of the desired connector in our 1D flattened array
        offset = sum([(self.connectors[c_idx].max_seq_length * 2 + 2) for i in range(0, c_idx)])

        # this connector also has a mu and sigma so we move 2 more indices forward
        offset += 2

        prob_sum = 0.0
        for i in range(len(p_list)):
            prob_sum += p_list[i]

            # here, the offset is used so that the correc indices in our 1D array of connector info
            # is modified

            # the line below corresponds the the pfs for the connector, and the line below it 
            # is for the aucs of the same connector
            c_scores[offset + i] = np.log(p_list[i])
            c_scores[offset + con.max_seq_length + i] = np.log(prob_sum)


    def get_placement(self, sequence: str) -> PlacementObject:
        """
        Calculates a placement for a given organism for a given sequence using
        the _multiplacement.calculate function.

        Preconditions:
            organism object has been flattend (see OrganismObject.flatten(self) function for details)

        Args:
            sequence: DNA sequence to place organism on

        Returns:
            PlacementObject containing information of optimal placement
        """
        # if the organism cannot be placed on the sequence, then it recieves a very low score
        # and is not attempted to be placed
        if self.sum_recognizer_lengths > len(sequence):
            raise ValueError(("This shouldn't have happened. The check_size " +
                              "function was supposed to prevent organisms from " +
                              "getting 'too big'. Something is not working."))
            # placement = Placement(self._id, sequence)
            # placement.set_energy(-1E100)
            # return placement
        
        num_recognizers = len(self.recognizers)
        
        # max length represents the gap length until precomputation has been done for, -1 means
        # no precomputation has been done
        max_length = -1
        
        # these arrays are allocated for the placement function to fill
        # gaps stores [starting position, gap_length_1, gap_length_2, ..., gap_length_n]
        # gap_scores stores the score associated with the gap length for the corresponding connector
        # recognizer_scores stores [rec_score_1, rec_score_2, ..., rec_score_n, total_score]
        # where the total score is the sum of the scores of all recognizers and connectors
        gaps = np.empty(num_recognizers, dtype = np.dtype('i'))
        gap_scores = np.empty(num_recognizers - 1, dtype = np.dtype('d'))
        recognizer_scores = np.empty(num_recognizers + 1, dtype = np.dtype('d'))
        
        # a temporary copy of the connector scores is made in case we need to rescale connector scores
        c_scores = np.copy(self.connectors_scores_flat)
        if num_recognizers > 1:
            
            #this determines until what gap length has been precomputed for [this will probably be removed since everything is precomputed now]
            max_length = self.connectors[0].max_seq_length - 1
            
            # this determines which connecters need to be rescaled and calls the function that rescales them
            # we rescale the scores for a connector when all pfs for gaps of length [0, sequence_length - minimum_length]
            # and mu > (length of sequence + 0.5)
            for i in range(len(self.connectors)):
                
                if float(len(sequence)) + 0.5 < self.connectors[i]._mu:
                    
                    # idx is the greatest possible gap that the connector could have in the placement
                    idx = len(sequence) - self.sum_recognizer_lengths
                    
                    # if the alternative model prediction for the gap of the maximum length is below our
                    # threshold then we rescale the scores for that connector
                    if self.connectors[i].stored_pdfs[idx] <= np.log(1E-15):

                        # the offset in the overall array of c scores will be the summation of the lengths of
                        # pfs, aucs, mu, and sigma for each connector
                        # each connector has max_seq_length pfs and aucs, and 1 for mu and 1 for sigma
                        # thats how we get (max_seq_length * 2 + 2) for each connector before the current one
                        # then 2 more indices for the current mu and sigma
                        offset = sum([(self.connectors[i].max_seq_length * 2 + 2) for i in range(0, i)]) + 2
                        self.connectors[i].adjust_scores(c_scores[offset:], len(sequence), self.sum_recognizer_lengths)
                        # XXX printing (for debugging purposes)
                        # print("rescaling...")
                        # print(len(c_scores))
                        
                        #self.adjust_connector_scores(i, c_scores, len(sequence))

        _multiplacement.calculate(bytes(sequence, "ASCII"), bytes(self.recognizer_types, "ASCII"), \
            self.recognizers_flat, self.recognizer_lengths,  c_scores, recognizer_scores, gap_scores, \
            gaps, max_length, self.recognizer_models, self.recognizer_bin_edges, self.recognizer_bin_nums)
        
        # parsing of the data from our placement function is done here
        # a placement object is instantiated and stores all of the scores that we obtained
        placement = PlacementObject(self._id, sequence)
        placement.set_energy(float(recognizer_scores[-1]))
        placement.set_recognizer_scores([float(score) for score in recognizer_scores[:-1]])
        placement.set_connectors_scores([float(score) for score in gap_scores])
        
        # the starting and ending positions for each recognizer and connector are calculated
        # using the starting position, lengths of each recognizer, and size of each gap
        current_position = gaps[0]
        placement.set_recognizer_types(self.recognizer_types)
        for i in range(num_recognizers):
            
            stop = current_position + self.recognizer_lengths[i]
            placement.append_recognizer_position([int(current_position), int(stop)])
            current_position += self.recognizer_lengths[i]
            
            if i < num_recognizers - 1:
                stop = current_position + gaps[i + 1]
                placement.append_connector_position([int(current_position), int(stop)])
                current_position += gaps[i + 1]
        
        return placement

    def set_sum_recognizer_lengths(self) -> None:
        self.sum_recognizer_lengths = sum([i.length for i in self.recognizers])
        
    def flatten(self) -> None:
        """
        Flattens an organisms recognizers into a 1D array
        computes gap scores and forms a 1D array out of them.
        This function is to be called whenever an organism is made
        or changed.


        Recognizers: 
            PSSMs are stored in a 1D array, self.recognizers_flat. The scores for
            each base are stored sequentially in the order [A, G, C, T] for each column
            of each PSSM.

            example layout:
        
               a PSSM with 3 columns

                 1    2    3
            A | x1A  x2A  x3A |
            G | x1G  x2G  x3G |  --> self.recognizers_flat = [x1A, x1G, x1C, x1T, x2A, x2G, x2C, x2T, x3A, x3G, x3C, x3T]
            C | x1C  x2C  x3C |
            T | x1T  x2T  x3T |

            Shape recognizers have a null model and an alternative model,
            both are stored in self.recognizer_models. This array holds these models
            for all shape recognizers in the organism, with the null model first, and
            alternative model second for each shape.

            example layout:

               arbitray Shape Recognizer with num_bins = n

            null_model: [N1, N2, N3, ... , Nn]
            alt_model:  [A1, A2, A3, ... , An]    -> self.recognizer_models = null_model + alt_model (simply the 2 arrays concatenated)
        
            bin_edges:  [E1, E2, E3, ... , En+1]  -> self.recognizer_bin_edges = bin_edges
            num_bins:   -> self.model_bin_nums = num_bins

            the size of each recognizer is stored in the array self.recognizer_lengths

            the type of a recognizer is represented by a char and all recognizer types are stored
            as a string, the char from each recognizer concatenated together, in self.recognizer_types

            p: PSSM - position specific scoring matrix
            m: mgw  - major groove width
            t: prot - propeller twist
            r: roll - roll
            h: helt - helical twist       
                
            Recognizers are arranged back to back in these arrays, directly one after another


        Connectors:
            a flattened connector is arranged into an 1D array storing its mu, sigma, pfs, and aucs.
            all of this information is held in self.connectors_scores_flat.

            pfs and aucs are computed for each interger value [0, max_seq_length] whenever mutated
            or created. The precomputed arrays store the log2 of the calculated probabilities

            example connector:
                _mu:    x
                _sigma: y
                stored_pdfs: [P0, P1, P2, ..., Pmax_seq_length]
                stored_cdfs: [C0, C1, C2, ..., Cmax_seq_length]

            Would look like [x, y, P0, P1, P2, ..., Pmax_seq_length, C0, C1, C2, ..., Cmax_seq_length]
            when flattend.

            Connectors are also stored back to back, so self.connectors_scores_flat would hold the information
            above for each connector in the whole organism.
            
        Args:
            None
        
        Returns:
            None
        """

        # instantiation of lists to hold flattened info
        flat_recognizers = []
        flat_connector_scores = []
        flat_rec_models = []
        rec_bin_nums = []
        rec_bin_edges = []
        recognizer_lengths = []
        self.recognizer_types = ""

        
        #each recognizer is iterated over and the required information is pulled and put into a 1D list
        for recognizer in self.recognizers:

            #if its a pssm we get the scores for each base for each column and put them into our
            #1D list
            if recognizer.get_type() == 'p':
                for column in recognizer.pssm:
                    for base in ['a','g','c','t']:
                        flat_recognizers.append(column[base])

            #if its a shape recognizer we get the information from its null_model, alt_model, and
            #put those into the same 1D list. And add the array containing the edges of each bin
            #to our edge array
            else:
                for prob in recognizer.null_model:
                    flat_rec_models.append(prob)
                for prob in recognizer.alt_model:
                    flat_rec_models.append(prob)
                for edge in recognizer.edges:
                    rec_bin_edges.append(edge)
                rec_bin_nums.append(len(recognizer.edges))
            recognizer_lengths.append(recognizer.length)
            self.recognizer_types += recognizer.get_type()

        # organism holds a numpy array of the flattened lists
        self.recognizers_flat = np.array(flat_recognizers, dtype = np.dtype('d')) 
        self.recognizer_lengths = np.array(recognizer_lengths, dtype = np.dtype('i'))
        self.recognizer_models = np.array(flat_rec_models, dtype = np.dtype('d'))
        self.recognizer_bin_nums = np.array(rec_bin_nums, dtype = np.dtype('i'))
        self.recognizer_bin_edges = np.array(rec_bin_edges, dtype = np.dtype('d'))

        # the minimum length of the organism is updated in case a recognizer was added, removed
        # or changed size
        self.set_sum_recognizer_lengths()


        if self.is_precomputed == True:

            #here we get the mu, sigma, pfs, and aucs for each connector and concatenate them all
            #into the same list
            for connector in self.connectors:
                flat_connector_scores += [connector._mu, connector._sigma]
                flat_connector_scores += (connector.stored_pdfs + connector.stored_cdfs)

            # organism holds a numpy array of the flattened lists
            self.connectors_scores_flat = np.array(flat_connector_scores, dtype = np.dtype('d'))

        # this is deprecated, precomputing is always better
        else:
            con_mu_sigma = []
            for connector in self.connectors:
                con_mu_sigma.append(connector._mu)
                con_mu_sigma.append(connector._sigma)

            self.connectors_scores_flat = np.array(con_mu_sigma, dtype=np.dtype('d'))
