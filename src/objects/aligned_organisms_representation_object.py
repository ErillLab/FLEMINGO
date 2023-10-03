# -*- coding: utf-8 -*-
"""
Aligned organisms representation object

"""

import copy

class AlignedOrganismsRepresentation:
    """
    Aligned organisms representation object
    
    """
    
    def __init__(self, organism1, organism2, dna_seq):
        """
        AlignedOrganismsRepresentation object constructor.

        Args:
            XXX
        """
        
        # Initialize organisms representations
        self.organism1 = []
        self.organism2 = []
        
        # Adjustment values for connectors
        self.org1_connectors_adjustments = None
        self.org2_connectors_adjustments = None
        
        # Independent units of recombination
        self.units = None
        
        # Connectors table (2D list). It stores at row i, column j the connector(s)
        # available to link index i to index j in the organism-representations
        self.connectors_table = None
        
        # Organisms IDs
        self.parents_ids = [organism1._id, organism2._id]
        self.children_ids = [None, None]
        
        # Compile attributes
        self.get_aligned_representations(organism1, organism2, dna_seq)
    
    def set_organsism1(self, nodes_list):
        self.organism1 = copy.deepcopy(nodes_list)
    
    def set_organsism2(self, nodes_list):
        self.organism2 = copy.deepcopy(nodes_list)
    
    def append_to_organism1(self, node_name):
        self.organism1.append(copy.deepcopy(node_name))
    
    def append_to_organism2(self, node_name):
        self.organism2.append(copy.deepcopy(node_name))
    
    def print_representation(self):
        print('\t'.join(self.organism1))
        print('\t'.join(self.organism2))
    
    def swap_unit(self, unit_start, unit_stop):
        '''
        Swaps the unit specified with the input bounds. The corresponding part
        of organism1 will end up into organism2, and vice versa.
        '''
        tmp = self.organism1[unit_start: unit_stop]
        self.organism1[unit_start: unit_stop] = self.organism2[unit_start: unit_stop]
        self.organism2[unit_start: unit_stop] = tmp
    
    def set_children_IDs(self, child1_id, child2_id):
        self.children_ids = [child1_id, child2_id]
    
    def get_aligned_representations(self, parent1, parent2, dna_seq):
        '''
        !!! Update docstring ...
        
        
        Places both the parents on the given DNA sequence, in order to 'align'
        them, one against the other. An abstract representation of the two
        aligned parents is returned as lists of symbols.
        
        EXAMPLE:
        This scheme
        
            p1_0    p1_1    -
            -       p2_0    p2_1
        
        says that recognizer 1 of parent1 ('p1') overlaps with recognizer 0 of
        parent2 ('p2'). Instead, recognizer 0 of parent 1 is unpaired, placing
        to the left of where parent2 is placed. Recognizer 1 of parent2 is also
        unpaired, placing to the right of where parent1 is placed.
        This scheme would be returned as a couple of lists:
        
            (
                ['p1_0', 'p1_1', '-'],
                ['-', 'p2_0', 'p2_1']
            )
        
        Returns
        -------
        parents_repres :
            Object of the class AlignedOrganismsRepresentation
        '''
        
        # Place the two organisms on the same DNA sequence
        placement1 = parent1.get_placement(dna_seq)
        placement2 = parent2.get_placement(dna_seq)
        
        # These dictionaries say which recognizer of an organism is occupying a
        # certain DNA position
        pos_to_recog_dict1 = self.get_pos_to_recog_idx_dict(placement1, 'p1')
        pos_to_recog_dict2 = self.get_pos_to_recog_idx_dict(placement2, 'p2')
        
        # Initialize symbolic representations
        p1_repres = []
        p2_repres = []
        
        # All the encountered pairs (each pair is made of one element from parent1,
        # and the other from parent2) are stored in this set
        pairs = set([])
        
        for i in range(len(dna_seq)):
            p1, p2 = '-', '-'
            
            if i in pos_to_recog_dict1.keys():
                p1 = pos_to_recog_dict1[i]
            
            if i in pos_to_recog_dict2.keys():
                p2 = pos_to_recog_dict2[i]
            
            pair = (p1, p2)
            # ignore DNA regions where there aren't recogs
            if pair != ('-','-'):
                
                # avoid repeating the match for all the DNA positions where the
                # match occurs
                if pair not in pairs:
                    pairs.add(pair)
                    # Compile parents representations
                    p1_repres.append(p1)
                    p2_repres.append(p2)
        
        # Remove protrusions
        p1_repres, p2_repres = self.remove_protrusions(p1_repres, p2_repres)
        
        # Set the attributes for organisms representations
        self.set_organsism1(p1_repres)
        self.set_organsism2(p2_repres)
        
        # Define units
        self.define_independent_units()
        
        # Set attributes: `org1_connectors_adjustments`, `org2_connectors_adjustments`
        self.set_connector_adjustments(placement1, placement2)
    
    
    def set_connector_adjustments(self, placement1, placement2):
        '''
        Sets two attributes:
            - self.org1_connectors_adjustments
            - self.org2_connectors_adjustments
        These two attributes are dictionaries where the keys are connector
        indexes (int) and the values are lists of two elements: the first is
        the adjustment for the left bound of the connector and the second is
        the adjustment for the right bound of the connector.
        
        EXAMPLE:
            {0: [0,0], 1: [0,0], 2:[-3, 1]}
            XXX ...
        '''
        
        # Initialize connector adjustments dictionaries
        org1_n_connectors = len(placement1.connectors_scores)
        org2_n_connectors = len(placement2.connectors_scores)
        org1_conn_adj_dict = {i : [0,0] for i in range(org1_n_connectors)}
        org2_conn_adj_dict = {i : [0,0] for i in range(org2_n_connectors)}
        
        # Compile connector adjustments dictionaries
        for unit_start, unit_stop in self.units:
            # Ignore if it's not an overlap between recognizers
            if (self.organism1[unit_start:unit_stop] == ["-"] or
                self.organism2[unit_start:unit_stop] == ["-"]):
                continue
            
            else:
                # LEFT connectors
                # ---------------
                p1recog_name = self.organism1[unit_start]
                p2recog_name = self.organism2[unit_start]
                p1recog_idx = int(p1recog_name.split('_')[1])
                p2recog_idx = int(p2recog_name.split('_')[1])
                # The indexes shouldn't be both 0 (there wouldn't be LEFT connectors)
                if p1recog_idx == 0 and p2recog_idx == 0:
                    pass
                else:
                    p1recog_start = placement1.recognizers_positions[p1recog_idx][0]
                    p2recog_start = placement2.recognizers_positions[p2recog_idx][0]
                    
                    # Store the right-bound adjustment for the connectors to the LEFT
                    displacement = p1recog_start - p2recog_start
                    
                    p1conn_idx = p1recog_idx - 1
                    if p1conn_idx >= 0:
                        org1_conn_adj_dict[p1conn_idx][1] -= displacement
                    
                    p2conn_idx = p2recog_idx - 1
                    if p2conn_idx >= 0:
                        org2_conn_adj_dict[p2conn_idx][1] += displacement
                
                # RIGHT connectors
                # ----------------
                p1recog_name = self.organism1[unit_stop - 1]
                p2recog_name = self.organism2[unit_stop - 1]
                p1recog_idx = int(p1recog_name.split('_')[1])
                p2recog_idx = int(p2recog_name.split('_')[1])
                # The indexes shouldn't be both the last (there wouldn't be RIGHT connectors)
                if p1recog_idx == org1_n_connectors and p2recog_idx == org2_n_connectors:
                    pass
                else:
                    p1recog_stop = placement1.recognizers_positions[p1recog_idx][1]
                    p2recog_stop = placement2.recognizers_positions[p2recog_idx][1]
                    
                    # Store the left-bound adjustment for the connectors to the RIGHT
                    displacement = p1recog_stop - p2recog_stop
                    
                    p1conn_idx = p1recog_idx
                    if p1conn_idx < org1_n_connectors:
                        org1_conn_adj_dict[p1conn_idx][0] += displacement
                    
                    p2conn_idx = p2recog_idx
                    if p2conn_idx < org2_n_connectors:
                        org2_conn_adj_dict[p2conn_idx][0] -= displacement
        
        # Set connector adjustments attributes
        self.org1_connectors_adjustments = org1_conn_adj_dict
        self.org2_connectors_adjustments = org2_conn_adj_dict    
    
    def remove_protrusions(self, p1_repres, p2_repres):
        '''
        Removes protrusions from the initial (draft) representations.
        
        Explanation:
        A 1-bp overlap between recognizers is enough for them to get 'paired'.
        This means that the overlap can be imperfect, with flanking parts of
        the recognizers being unpaired. Those will be ignored.
        
        EXAMPLE:
        If two recognizers are placed on DNA in this way
            -----AAAA-------
            -------BBBB-----
        the desired representation is
            A
            B
        and not
            AA-
            -BB
        Therefore, the two positions to the left and to the right of the
        A-B match will be called 'protrusions', and they will be removed.
        '''
        matches_p1 = set([])  # recogs of parent1 that overlap with a recog
        matches_p2 = set([])  # recogs of parent2 that overlap with a recog
        
        for i in range(len(p1_repres)):
            p1_node, p2_node = p1_repres[i], p2_repres[i]
            if p1_node != '-' and p2_node != '-':
                matches_p1.add(p1_node)
                matches_p2.add(p2_node)
        
        # If recognizer X is in matches_p1 or matches_p2 (menaing that it
        # overlaps at least once with another recognizer) all the other
        # eventual pairings of X with "-" are protrusions.
        
        # Here we store the indexes of the positions where there's a 'protrusion'
        protrusions = []        
        
        for i in range(len(p1_repres)):
            if p1_repres[i] in matches_p1:
                if p2_repres[i] == '-':
                    protrusions.append(i)
        
        for i in range(len(p2_repres)):
            if p2_repres[i] in matches_p2:
                if p1_repres[i] == '-':
                    protrusions.append(i)
        
        # Skip the protrusions and return the desired representations
        p1_repres = [p1_repres[i] for i in range(len(p1_repres)) if i not in protrusions]
        p2_repres = [p2_repres[i] for i in range(len(p2_repres)) if i not in protrusions]
        
        return p1_repres, p2_repres
    
    def get_pos_to_recog_idx_dict(self, placement, org_tag):
        '''
        For the given `placement`, each DNA position covered by some recognizer
        is mapped to a string that says what recognizer is placed there. The
        string will contain a tag for the organism (`org_tag` argument), joint
        with the recog index by an underscore.
        
        EXAMPLE:
        This dictionary
            {112: 'p1_0', 113: 'p1_0', 114: 'p1_0', 115: 'p1_0'}
        will be used to know that (for example) position 114 is covered by
        recognizer 0. In this case, 'p1' was the org_tag value specified as
        input, used to identify an organism.
        
        Parameters
        ----------
        placement : PlacementObject
        org_tag : string

        Returns
        -------
        pos_to_recog_dict : dictionary
        '''
        recog_positions = placement.recognizers_positions
        
        pos_to_recog_dict = {}
        
        # for each recognizer
        for i in range(len(recog_positions)):
            # start and stop DNA positions of recognizer i
            start, stop = recog_positions[i]
            
            # DNA positions occupied by this recognizer
            for pos in range(start, stop):
                pos_to_recog_dict[pos] = org_tag + '_' + str(i)  # i is the recog idx
        
        return pos_to_recog_dict
    
    
    def annotate_available_connectors(self):
        '''
        Sets the `connectors_table` attribute.
        The representations of the aligned parents are lists of symbols, stored
        in the `organism1` and `organism2` attributes.
        For each possible couple of positions in the representations, this
        function annotates whether the parents have a connector that connects
        them. This info is stored in the attribute `connectors_table`.
        
        EXAMPLE:
        In this representations
        
            p1_0    p1_1    -       p1_2
            -       p2_0    p2_1    -   
        
        index 0 is liked to index 1 by the first connector of p1 (parent1): the
        connector that connects recognizer p1_0 with recognizer p1_1.
        
        Index 1 is liked to index 3 by the second connector of p1: the
        connector that connects recognizer p1_1 with recognizer p1_2.
        
        Index 1 is liked to index 2 by the only connector of p2: the
        connector that connects recognizer p2_0 with recognizer p2_1.
        
        Sets `connectors_table`: a table in the form of a 2D list.
        This table stores at row i, column j the connector(s) available to link
        index i to index j.
        '''
        
        # !!! Add code here to catch case when organisms were not set yet
        n = len(self.organism1)
        
        # 2D list where each item is an emtpy list
        connectors_table = [[ [] for i in range(n)] for j in range(n)]
        
        # Each parent representation is coupled with a tag ('p1' or 'p2')
        parents = [(self.organism1, 'p1'), (self.organism2, 'p2')]
        
        for (org_repr, org_tag) in parents:
            
            # Indexes where a recognizer of this parent is present        
            recogs_indexes = []
            for idx in range(len(org_repr)):
                if org_repr[idx] != '-':
                    recogs_indexes.append(idx)
            
            connector_idx = 0
            for i in range(len(recogs_indexes)-1):
                
                left_recog_idx = recogs_indexes[i]
                right_recog_idx = recogs_indexes[i+1]
                
                left_recog_name = org_repr[left_recog_idx]
                right_recog_name = org_repr[right_recog_idx]
                
                if left_recog_name != right_recog_name:                
                    connector_name = org_tag + '_' + str(connector_idx)
                    connector_idx += 1
                    connectors_table[left_recog_idx][right_recog_idx].append(connector_name)
        
        # Set the table storing info about available connectors
        self.connectors_table = connectors_table
    
    def define_independent_units(self):
        '''
        This function is used to define what chunks of the organisms'
        representation are going to work as independent units in the
        recombination process. Within each unit, the part from parent1 can be
        swapped with the part from parent2 (by get_aligned_children_repr function)
        with 50% probability.
        This function generates a list of units' spans, where each element is a
        (start, stop) tuple.
        
        EXAMPLE
        In this representations
        
            p1_0    p1_1    -       p1_2
            p2_0    p2_0    p2_1    -   
        
        p2_0 is partially overlapping with p1_0, and partially with p1_1. In
        this scenario, the first two positions of the representations work as a
        single unit. Therefore, the obtained list of independent units will be
        
            [ (0, 2), (2, 3), (3, 4) ]
        
        '''
        
        # !!! Add code here to catch case when organisms were not set yet
        org1_repr = self.organism1
        org2_repr = self.organism2
        
        # Initialize list about where each unit starts
        unit_starts = [0]
        
        for i in range(1, len(org1_repr)):
            
            # If in org1_repr at position i there is the same recognizer as the
            # one at position i-1
            if org1_repr[i] != '-' and org1_repr[i] == org1_repr[i-1]:
                # Then this is not yet the start of the next unit
                continue
        
            # If in org2_repr at position i there is the same recognizer as the
            # one at position i-1
            if org2_repr[i] != '-' and org2_repr[i] == org2_repr[i-1]:
                # Then this is not yet the start of the next unit
                continue
            
            unit_starts.append(i)  # i is the start position of a new unit
        
        # Each unit stops where the next unit starts (or when the list ends in
        # the case of the last unit)
        unit_stops = unit_starts[1:] + [len(org1_repr)]
        
        # Make a list of units. Each unit is a tuple: (start, stop)
        units = list(zip(unit_starts, unit_stops))
        self.units = units










