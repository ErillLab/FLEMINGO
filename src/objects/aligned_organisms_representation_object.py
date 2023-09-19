# -*- coding: utf-8 -*-
"""
Aligned organisms representation object

"""

import copy

class AlignedOrganismsRepresentation:
    """
    Aligned organisms representation object
    
    """
    
    def __init__(self, organism1_id, organism2_id):
        """
        AlignedOrganismsRepresentation object constructor.

        Args:
            XXX
        """
        
        # Initialize organisms representations
        self.organism1 = []
        self.organism2 = []
        # XXX
        self.org1_connectors_adjustments = None
        self.org2_connectors_adjustments = None
        
        # XXX
        self.units = None
        
        # XXX
        self.connectors_table = None
        
        # Organisms IDs
        self.parents_ids = [organism1_id, organism2_id]
        self.children_ids = [None, None]
    
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
        # XXX
        # Perform the swapping, which means that the part from organism1 will
        # end up into organism2, and the part from organism2 will end up into
        # organism1.
        tmp = self.organism1[unit_start: unit_stop]
        self.organism1[unit_start: unit_stop] = self.organism2[unit_start: unit_stop]
        self.organism2[unit_start: unit_stop] = tmp
    
    def set_children_IDs(self, child1_id, child2_id):
        self.children_ids = [child1_id, child2_id]
    
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
        
        # XXX Code here to catch case when organisms were not set yet
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
        
        # XXX Code here to catch case when organisms were not set yet
        org1_repr = self.organism1
        org2_repr = self.organism2
        
        # Initialize list about where each unit starts
        unit_starts = [0]
        
        for i in range(1, len(org1_repr)):
            
            # If in org1_repr at position i there is the same recognizer as the one
            # at position i-1
            if org1_repr[i] != '-' and org1_repr[i] == org1_repr[i-1]:
                # Then this is not yet the start of the next unit
                continue
        
            # If in org2_repr at position i there is the same recognizer as the one
            # at position i-1
            if org2_repr[i] != '-' and org2_repr[i] == org2_repr[i-1]:
                # Then this is not yet the start of the next unit
                continue
            
            unit_starts.append(i)  # i is the start position of a new unit
        
        # Each unit stops where the next unit starts (or when the list ends in the
        # case of the last unit)
        unit_stops = unit_starts[1:] + [len(org1_repr)]
        
        # Make a list of units. Each unit is a tuple: (start, stop)
        units = list(zip(unit_starts, unit_stops))
        self.units = units









