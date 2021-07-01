# -*- coding: utf-8 -*-
"""
Organism object
It allocates the full data structure
"""


import random
import numpy as np
from scipy.stats import ks_2samp
import copy
from .placement_object import PlacementObject


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
	
        # assign organism-specific parameters
	
        # whether fitness is computed over sequences as sum or average
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
        
    def mutate(self, org_factory) -> None:
        """Mutates an organism based on JSON configured probabilities

        Args:
            org_factory (OrganismFactory): Factory of organisms and node
                                           components
        """

        
        # Delete a recognizer (and one parent connector)
        if random.random() < self.mutate_probability_delete_recognizer:
            
            n_recognizers = self.count_recognizers()
            # if this is a single-node organism, skip deletion
            if n_recognizers != 1:
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
                    
    				# recreate vector of recognizers and connectors
    				# skipping the recgonizer/connector selected for deletion
                    new_recognizers = (self.recognizers[:recognizer_idx] +
                                       self.recognizers[recognizer_idx + 1:])
                    new_connectors = (self.connectors[:connector_idx] +
                                       self.connectors[connector_idx + 1:])
                    
    				# assign new vectors
                    self.recognizers = new_recognizers
                    self.connectors = new_connectors
                
    			# "intelligent" method: a new connector is created merging
    			# the information of the right and left connectors for the
    			# recognizer targeted for deletion
                if self.deletion_method == "intelligent":
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
    
    				# recreate vector of recognizers and connectors
    				# skipping the recgonizer/connector selected for deletion					
                    new_recognizers = (self.recognizers[:recognizer_idx] +
                                       self.recognizers[recognizer_idx + 1:])
                    new_connectors = (self.connectors[:connector_idx] +
                                       self.connectors[connector_idx + 1:])
                    
    				# assign new vectors
                    self.recognizers = new_recognizers
                    self.connectors = new_connectors
        
        
        # Insert a recognizer (and one parent connector)
        if random.random() < self.mutate_probability_insert_recognizer:
            
			# instantiate the new recognizer and connector
            new_connector = org_factory.create_connector()
            new_recognizer = org_factory.create_pssm()
            
			# "blind" method: one of the existing connectors is used
			# (unchanged) to connect to the new recognizer
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
                
				# recreate vector of connectors/recognizers, adding
				# the newly minted recognizer+connector
                new_recognizers = (self.recognizers[:recognizer_idx] +
                                   [new_recognizer] +
                                   self.recognizers[recognizer_idx:])
                new_connectors = (self.connectors[:connector_idx] +
                                   [new_connector] +
                                   self.connectors[connector_idx:])
                
				# assign new connector/recognizer vectors
                self.recognizers = new_recognizers
                self.connectors = new_connectors
            
            
			# "intelligent" method: the existing and new connector are
			# "averaged" so that their mu and sigma are "equivalent" to
			# to the connector previously occupying the position where
			# the new recognizer has been inserted
            if self.insertion_method == "intelligent":
                
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

					# if the new recognizer is NOT be the first in the chain
                    if recognizer_idx != 0:
                        connector_to_compress = recognizer_idx - 1
                        # (No compression is required if insertion occurs to
                        # the left of the first recognizer of the chain)
                        
                else:  # Insertion occurs to the right
                    # First connector after insertion point
                    connector_idx = recognizer_idx
                    
					# if the new recognizer is NOT be the last in the chain
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
                    by the old recognizer, so that the insertion doesn't heavily
                    affect the placement of the nodes to the sides of the
                    insertion point. So we need to scale down the mus and sigmas
                    accordingly.'''
                    
                    # Adjust MUs
                    ''' If the legth of the inserted recognizer alone is larger
                    than the gap (i.e. larger than the mu of the already present
                    connector) we can't reach the goal entirely and the best
                    thing we can do is to keep the mus of the two connectors
                    maximally shrinked. Otherwise, the mus of the connectors
                    will be scaled so that their sum is equal to the gap (mu of
                    the old connector) minus the length of the inserted
                    recognizer.'''
                    
                    mu_old_connector = self.connectors[connector_to_compress]._mu
                    mu_new_connector = new_connector._mu                    
                    current_sum_mus = mu_old_connector + mu_new_connector
                    expected_sum_mus = mu_old_connector - new_recognizer.length
                    # If the inserted recognizer alone is larger than the gap
                    if expected_sum_mus < 0:
                        # best thing we can do is to maximally shrink mus
                        expected_sum_mus = 0
                    
                    if current_sum_mus == 0:  # To avoid 0/0 case
                        mu_scaling_factor = 1
                    else:
                        mu_scaling_factor =  expected_sum_mus / current_sum_mus
                    # Compress the neighbour connector
                    self.connectors[connector_to_compress].set_mu(
                        mu_old_connector * mu_scaling_factor
                    )
                    # Compress the new inserted connector
                    new_connector.set_mu(
                        mu_new_connector * mu_scaling_factor
                    )
                    
                    # Adjust SIGMAs
                    ''' Variances are additive (under assumption of independence
                    between the two random variables). Therefore, the overall
                    variance will be the sum of the variances of the two
                    connectors. The variances will be scaled so that their sum
                    is equal to the variance of the connector that was already
                    present before the insertion. The standard deviations,
                    which are the square root of the variances, will be adjusted
                    accordingly.'''
                    var_old_connector = self.connectors[connector_to_compress]._sigma ** 2
                    var_new_connector = new_connector._sigma**2
                    current_sum_variances = var_old_connector + var_new_connector
                    expected_sum_variances = var_old_connector
                    
                    if current_sum_variances == 0:  # To avoid 0/0 case
                        var_scaling_factor = 1
                    else:
                        var_scaling_factor =  expected_sum_variances / current_sum_variances
                    # Compress the neighbour connector
                    self.connectors[connector_to_compress].set_sigma(
                        np.sqrt(var_old_connector * var_scaling_factor)
                    )
                    # Compress the new inserted connector
                    new_connector.set_sigma(
                        np.sqrt(var_new_connector * var_scaling_factor)
                    )
                
				# recreate vector of connectors/recognizers, adding
				# the newly minted recognizer+connector and containing
				# also the "compressed" existing connector
                new_recognizers = (self.recognizers[:recognizer_idx] +
                                   [new_recognizer] +
                                   self.recognizers[recognizer_idx:])
                new_connectors = (self.connectors[:connector_idx] +
                                   [new_connector] +
                                   self.connectors[connector_idx:])
                
				# assign new connector/recognizer vectors
                self.recognizers = new_recognizers
                self.connectors = new_connectors
        
        
        # Mutate nodes
        # If MUTATE_PROBABILITY_NODE_MUTATION is set to real value, a single
        # node is selected. If set to null, all nodes are mutated
        if self.mutate_probability_node_mutation:
            # Mutate a random node
            if random.random() < self.mutate_probability_node_mutation:
    
                n_nodes = self.count_nodes()
                random_node_idx = random.randint(0, n_nodes - 1)
                if random_node_idx < self.count_recognizers():
                    # mutate a recognizer
                    self.recognizers[random_node_idx].mutate(org_factory)
                else:
                    # mutate a connector
                    connector_idx = random_node_idx - self.count_recognizers()
                    self.connectors[connector_idx].mutate(org_factory)
        else:
            for recognizer in self.recognizers:
                recognizer.mutate(org_factory)
            for connector in self.connectors:
                connector.mutate(org_factory)
        
        # no matter what mutation is applied, order/columns/number of pssm's
        # may have changed, so we call the set_row_to_pssm to set their mapping
        # on the alignment matrix anew
        self.set_row_to_pssm()
    
    def get_placement(self, dna_sequence, traceback=False) -> dict:
        """Places the organism elements (recognizers and connectors) on a sequence
		   in an optimal way, maximizing the energy (i.e. cumulative scores) obtained.
		   
		   That is, it returns the best possible placement of the organism on the
		   sequence, as a function of the cumulative organism energy.
		   
		   Inputs:
		   - dna_sequence: DNA sequence to place on
		   - print_out: bool to indicate whether to print or not the placement
           - out_file: file handle to write placement to, if desired
		
		   The placement function implements a modified Needleman-Wünch algorithm.
		   
		   The algorithm uses a two-dimensional matrix, with the sequence on the X-axis
		   and the PSSM columns on the Y-axis.
		   - Substitution scores are computed as PSSM column scores.
		   - First row is set to zeros
		   - First column is set to -inf
		   - Gaps are only allowed in terminal PSSM rows (last column of PSSM)
		     - Gaps can only take place between end of PSSM (cell) and another cell in row
		   - Gaps are scored according to the score from the connector between the two
		     PSSMs, provided with the distance (and its internal mu and sigma)
		   - Contiguous (diagonal) PSSM alignment invokes a zero gap using the appropriate
		     connector
		   - The alignment matrix is (M+1)x(N+1)
		   - The traceback matrix (which stores info for traceback) is 2 x (M+1)x(N+1).
		     - The extra dimension captures the two coordinates for the traceback to cell
		   - Traceback through a gap enforces that a diagonal move must be taken next
		     (this avoids double gaps in a row)
		   - Traceback is initiated at the cell with the best value on the bottom row
		"""
    
        # Initialize the two matrices (alignment + traceback matrices)
        
        # Number of rows
        m = self.sum_pssm_lengths()
        # Number of columns
        n = len(dna_sequence)
        
        # Initialize matrix of scores (alignment matrix)
		# Matrix is (M+1)x(N+1), with the extra "fake" row/columns
		# Matrix is initialized to -inf, then first row set to zero
        scores_matrix = np.full((m+1, n+1), -1 * np.inf)
        scores_matrix[0,:] = 0
        
        # Initialize matrix of pointers (traceback matrix), to None
		# Matrix is 2x(M+1)x(N+1), to store row/col of incoming cell
        pointers_matrix = np.full((2, m+1, n+1), None)
        
        # Fill the matrices (top-to-bottom, then left-to-right)
        for i in range(1, m + 1):
			# Row fill up is done in two passes:
			#  - First fill up with possible diagonal scores
			#  & (for terminal recognizer rows only)
			#  - Fill up with possible gap scores (horizontal moves)
			
            # Diagonal scores over row i
            for j in range(1, n + 1):
                # call PSSM column score function with the DNA sequence
                diag_score = self.get_diag_score(pointers_matrix, i, j, dna_sequence)
				# assign cumulative score to alignment matrix
                scores_matrix[i,j] = scores_matrix[i-1, j-1] + diag_score
                # Annotate "where we came from" in the pointers_matrix
                pointers_matrix[0][i,j] = i - 1  # row idx of the origin
                pointers_matrix[1][i,j] = j - 1  # column idx of the origin
                
            # Horizontal scores over row i
            # (only in rows at the interface with the next PSSM)
            if self.is_last(i) and i != m:
                
                # Scores and pointers from horizontal moves are temporarily stored in
                # the following arrays. They will be written altogether at the end of
                # the loop over j, to avoid reading them as starting scores for other
                # horizontal moves at a later cycles in the for loop over j
				
				# That is, the matrix with diagonal scores is left untouched, and used
				# as reference for any gap evaluations. The result of such gap evaluations
				# is placed on a temporary matrix. The best gap scores are computed there.
				# This is to avoid adding the gap-to-gap score (instead of
				# the diagonal score) that has replaced a diagonal score, as we move
				# further right (because otherwise you could "carry" multiple gap
				# scores
                tmp_gap_scores = scores_matrix[i,:].copy()  # vector of length n+1
                tmp_gap_pointers = pointers_matrix[:,i,:].copy()  # 2 x (n+1) matrix
                
				# for each column of the matrix
                for j in range(1, n + 1):
                    # Compute all the possible values from all the possible
                    # horizontal moves (gaps) that land on [i,j]
                    for start in range(j):
                        gap_size = j - start	# obtain distance for connector
						# obtain PSSM index, which is also the connector index
						# because this is the preceding PSSM
                        pssm_idx = self.row_to_pssm[i][0]
						# evaluate the score for connector given distance
						# the DNA sequence length (n) is needed for the connector
						# background model
                        gap_score = self.get_gap_score(pssm_idx, gap_size, n)
						# add cumulative score (using the score matrix, which has
						# not been modified by horizontal scores)
                        candidate_score = scores_matrix[i, start] + gap_score
                        
						# assess whether temp scores should be overwritten
						# that is, whether this gap is better than other gaps
                        if candidate_score >= tmp_gap_scores[j]:
                            # Write horizontal score if better than existing
                            tmp_gap_scores[j] = candidate_score
                            # Annotate "where we came from" in the tmp_gap_pointers
                            tmp_gap_pointers[0,j] = i  # row idx of the origin
                            tmp_gap_pointers[1,j] = start  # column idx of the origin
                
                
                # Update the original matrices
				# tmp_gap_scores contains the updated matrix, with the diagnonal moves
				# and any best horizontal moves replacing them if adequate
                scores_matrix[i,:] = tmp_gap_scores
                pointers_matrix[:,i,:] = tmp_gap_pointers
            
        # Get best binding energy (max value on bottom row)
        last_row = scores_matrix[-1,:]
        best = max(last_row)
        
        # BACKTRACKING
        
        # Initialize placement object
        placement = PlacementObject(self._id, dna_sequence)
        
        # Set the total binding energy in the placement object
        # Applying lower bound to energy if required
        if self.energy_threshold_method == "organism":
            E_threshold_value = self.energy_threshold_value
            if best < E_threshold_value:
                placement.set_energy(E_threshold_value)
            else:
                placement.set_energy(best)
        
        
        if traceback:
            # Position of best (where backtracking starts from)
            best_i = m  # it always comes from the last row by definition
            # if multiple positions in last row have best value, pick first
            # to initiate traceback
            best_j = int(np.where(last_row == best)[0][0])  # column of best value
            
            # Traverse back the matrix from the best element in the last row
            # and store the alignment path
    		# traverse_matrix is a recursive function that will generate the path
    		# taken by the optimal alignment
            alignment_path = []
            alignment_path = self.traverse_matrix(pointers_matrix, best_i, best_j, alignment_path)
            alignment_path.reverse()  # Top-down instead of bottom-up
            
            # Get scores and positions of all the nodes of the organism
            node_scores, node_pos_ends, cols_of_0_gaps = self.get_node_positions_and_energies(
                alignment_path, scores_matrix, pointers_matrix, dna_sequence
            )

            # Split node-scores into recognizers-scores and connectors-scores
            # Remove token node [first row is treated as a node, to provide
            # a start position for the following recognizer, by tracking its
            # end]
            node_scores = node_scores[1:]  
            recognizers_scores = []
            connectors_scores = []
            for i in range(len(node_scores)):
                if i % 2 == 0:
                    recognizers_scores.append(node_scores[i])
                else:
                    connectors_scores.append(node_scores[i])
            
            # Set the individual node energies in the placement object
            placement.set_recognizer_scores(recognizers_scores)
            placement.set_connectors_scores(connectors_scores)
            
            # Annotate where each node is placed (start and stop position on DNA)
            # in the placement object
            
            for i in range(0, len(node_pos_ends)-1):  # even indexes
                
                # get sequence coordinates for the recognizer
                start = node_pos_ends[i]
                stop = node_pos_ends[i+1]
                
                # Detect a post 0-bp gap recognizer (special case)
                if start in cols_of_0_gaps:
                    start -= 1
                
                # Define position of the node on DNA as a tuple
                node_position = (start, stop)
                
                if i % 2 == 0:
                    placement.append_recognizer_position(node_position)
                else:
                    placement.append_connector_position(node_position)
        
        return placement
    
    def get_additive_fitness(self, a_dna: list) -> dict:
        """

        Args:
            a_dna: list of dna sequences

        Returns:
            average/sum of the energy of the organism on the sequences
        """

        scores = []
		# for each sequence in the provided sequence set
        for s_dna in a_dna:
            placement = self.get_placement(s_dna)
            energy = placement.energy  # energy          
            scores.append(energy)
        
        score_stdev = np.std(scores)
        if self.cumulative_fit_method == "sum":
            # Compute fitness score as sum over the positive scores
            score = np.sum(scores)
        
        elif self.cumulative_fit_method == "mean":
            # Compute fitness score as average positive score
            score = np.mean(scores)
            
        elif self.cumulative_fit_method == "median":
            # Compute fitness score as median positive score
            score = np.median(scores)
        
        return {"score": score, "stdev" : score_stdev}
    
    def get_binding_energies(self, a_dna: list, traceback=False) -> list:
        """Return the binding energies for an array of DNA sequences.

        Args:
            a_dna: list of dna sequences

        Returns:
            list of binding eneregies
        """

        binding_energies = []
		# for each sequence in the provided sequence set
        for s_dna in a_dna:
            placement = self.get_placement(s_dna)
            energy = placement.energy
            binding_energies.append(energy)
        
        return binding_energies

    def get_kolmogorov_fitness(self, pos_dataset: list, neg_dataset: list,
                               traceback=False) -> float:
        """Returns the organism's fitness, defined as the Kolmogorov-Smirnov
           test statistic. This is bounded in [0,1].
           Test null assumes the samples are drawn from the same (continuous)
           distribution.
           The statistic is sensitive to differences in both location and shape 
           of the empirical cumulative distribution functions of the two samples.
        Args:
            pos_dataset: list of dna sequences in the positive dataset
            neg_dataset: list of dna sequences in the negative dataset
        Returns:
            fitness assigned to the organism
        """       
        # Values on the positive set
        pos_values = []
        for s_dna in pos_dataset:
            placement = self.get_placement(s_dna)
            pos_values.append(placement.energy)  # append sequence energy score
        
        # Values on the negative set
        neg_values = []
        for s_dna in neg_dataset:
            placement = self.get_placement(s_dna)
            neg_values.append(placement.energy)  # get sequence energy score
        
        # Compute fitness score as a Boltzmannian probability
        kolmogorov_fitness = ks_2samp(pos_values, neg_values).statistic
        
        return {"score": kolmogorov_fitness}
    
    def get_boltz_fitness(self, pos_dataset: list, neg_dataset: list,
                          genome_length: int) -> float:
        """Returns the organism's fitness, defined as the probability that the regulator binds a
        positive sequence. All the binding energies are turned into probabilities according to a
        Boltzmannian distribution. The probability of binding a particular sequence, given the binding
        energy on that sequence, is p = e**binding_energy / Z
        where Z is the partition function.
        A high number of negative sequences is assumed to be present (emulating the environment of a
        regulator that needs to find its targets on an entire genome).
        A coefficient called neg_factor is computed, so that the value of Z can be as high as if there
        were as	many negative sequences as required to cover the entire genome.

        Args:
            pos_dataset: list of dna sequences in the positive dataset
            neg_dataset: list of dna sequences in the negative dataset
            genome_length: integer representing the length of the genome

        Returns:
            fitness assigned to the organism
        """
        
        # Values on the positive set
        pos_values = []
        for s_dna in pos_dataset:
            placement = self.get_placement(s_dna)
            boltz_exp = np.e**placement.energy  # exp(energy)
            pos_values.append(boltz_exp)
        
        # Values on the negative set
        neg_values = []
        neg_lengths = []
        for s_dna in neg_dataset:
            placement = self.get_placement(s_dna)
            boltz_exp = np.e**placement.energy  # exp(energy)
            neg_values.append(boltz_exp)
            neg_lengths.append(len(s_dna))
        
        # Scaling factor, used to over-represent the negative scores, so that
        # it simulates a genome of specified length
        neg_factor = genome_length//sum(neg_lengths)
        
        # Partition function
        Z = sum(pos_values) + neg_factor * sum(neg_values)
        
        # Compute fitness score as a Boltzmannian probability
        boltz_fitness = sum(pos_values) / Z
        
        return {"score": boltz_fitness}

    def count_nodes(self) -> int:
        """Returns the number of nodes of the organism

        Returns:
            Number of nodes of the organism
        """

        return 2 * len(self.recognizers) - 1
    
    def count_connectors(self) -> int:
        """Returns the number of connectors of the organism

        Returns:
            Number of connectors.
        """
        
        return len(self.connectors)

    def count_recognizers(self) -> int:
        """Returns the number of recognizers of the organism

        Returns:
            Number of recognizers.
        """
        
        return len(self.recognizers)
    
    def sum_pssm_lengths(self) -> int:
        """Returns the sum of the lengths of all the PSSMs of the organism.
        """
        
        sum_lengths = 0
        for pssm in self.recognizers:
            sum_lengths += pssm.length
        
        return sum_lengths
    
    def get_gap_score(self, connector_idx, d, s_dna_len):
        """Calls the appropriate connector, with the given distance and length of the DNA sequence to
		   obtain the energy of the connector.
		"""
        if d == s_dna_len:
            return -1 * np.inf
        
        gap_score = self.connectors[connector_idx].get_score(d, s_dna_len)
        return gap_score
    
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
            zero_gap_score = connector.get_score(0, len(dna_sequence))
            diag_score += zero_gap_score
        
		# get nucleotide and compute PSSM score for it
        nucleotide = dna_sequence[col_idx - 1] 
        pssm_score = self.get_score_from_pssm(row_idx, nucleotide)
        diag_score += pssm_score
        
        return diag_score
    
    def is_first(self, row_idx_from_placement_matrix):
        """Returns true if we are on the first element of a PSSM recognizer
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
    
    def traverse_matrix(self, pointers_mat, i, j,
                        alignment_path=[], from_gap_flag=False) -> list:
        """Recursive function used for traceback
		   - i and j are the starting positions
		   - alignment_path is the path that is filled up recursively
		   - from_gap_flag identifies the case that we got to this position from a gap
		     and we therefore MUST go diagonal (no gap concatenation is allowed)
		"""
        
        # End of the recursion
        if i == None:  # the top row has been reached
            return alignment_path
        
		# add current cell coordinates to path
        alignment_path.append([i,j])
        
		# avoid double gaps
        if from_gap_flag == True:
            # move diagonally
            i_next = i - 1
            j_next = j - 1
	    # or move back through annotated pointer
        else:
            # move where the pointers-matrix says
            i_next = pointers_mat[0, i, j]
            j_next = pointers_mat[1, i, j]
        
		# if moving horizontally, indicate that with flag
        if i_next == i:
            return self.traverse_matrix(pointers_mat, i_next, j_next,
                                        alignment_path, from_gap_flag=True)
        else:
            return self.traverse_matrix(pointers_mat, i_next, j_next,
                                        alignment_path,from_gap_flag=False)

    
    def get_node_positions_and_energies(self, alignment_path, scores_matrix,
                                        pointers_matrix, dna_seq) -> list:
        """Takes the alignment path and the completed score and pointer matrices.
           Returns a list that contains:
           - node_scores: score of recognizer/connector node
           - node_placements_right_ends: column of matrix where placement of node ends
           - columns_of_0_bp_gaps: column with special contiguous recognizer case
           
           The alignment path is already reversed, so first element is top-left.
           
           Function goes through the alignment path and reports where each node
           ends (column in alingment matrix)
		"""
        
        # Use the path to get info about individual positions and score of all the
        # nodes of the organism
        
        node_scores = []
        previous_score = 0
        node_placements_right_ends = []
        previous_element = [None, None]
        columns_of_0_bp_gaps = []
        
        # for each element of the alignment path
        for element in alignment_path:
            row, column = element
            
            # if we are on a row where gaps are allowed (PSSM ends or connector)
            if self.is_last(row):
                
                # if we just landed on this row via a diagonal move, then 
                # this is a PSSM end
                if previous_element[0] == row - 1:
                    # The PSSM score needs to be recomputed (it could have been
                    # over-written by a gap-score)
                    
                    score = self.get_diag_score(pointers_matrix, row,
                                                         column, dna_seq)
                    cumulative_score = scores_matrix[row-1, column-1] + score
                
                # this was a gap, so no need to recompute score
                else:
                    cumulative_score = scores_matrix[row, column]
                
                # compute the node score, by substracting from cumulative
                node_score = cumulative_score - previous_score
                node_scores.append(node_score)
                previous_score = cumulative_score
                
                # mark its last position on matrix (right end)
                node_placements_right_ends.append(column)
                
            # if we are on the first position of a new pssm which is adjacent to
            # the previous one (gap of 0 bp), the cumulative score we read also
            # contains the score of the first position of the PSSM (not only the
            # connector score)
            if self.is_a_0_bp_gap(pointers_matrix, row, column):
                
                
                nucleotide = dna_seq[column - 1] 
                pssm_contribution = self.get_score_from_pssm(row, nucleotide)
                
                
                cell_score = scores_matrix[row, column]
                # Remove the pssm contribution, so that the cumulative score
                # doesn't include the first position of the next PSSM, but only the
                # contribution from the connector
                cumulative_score = cell_score - pssm_contribution
                node_score = cumulative_score - previous_score
                node_scores.append(node_score)
                previous_score = cumulative_score
                
                node_placements_right_ends.append(column)
                
                columns_of_0_bp_gaps.append(column)
            
            
            previous_element = [row, column]
            
        
        return [node_scores, node_placements_right_ends, columns_of_0_bp_gaps]
    
    
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
        organism_file = open(filename, "w+")
        organism_file.write("***** Organism {} *****".format(self._id))
        
        for i in range(len(self.recognizers) - 1):
            self.recognizers[i].export(organism_file)
            self.connectors[i].export(organism_file)
        self.recognizers[-1].export(organism_file)

        organism_file.write("\n")
        organism_file.close()

    def export_results(self, a_dna: list, filename: str) -> None:
        """Exports the binding profile of the organism against each of the 
           DNA sequences provided as a list

        Args:
            filename: Name of the file to export sequences
            a_dna: list fo sequences to export

        """
        
        ofile = open(filename, "w")
        # for each DNA sequence
        for s_dna in a_dna:
            placement = self.get_placement(s_dna.lower(), traceback=True)
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
        placement = self.get_placement(s_dna.lower(), traceback=True)
        placement.print_placement(stdout = True)

