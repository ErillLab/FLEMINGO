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
            
        """
        
        # Initialize organisms representations
        self.organism1 = []
        self.organism2 = []
        
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
        The representations of the aligned parents are lists of symbols.
        For each possible couple of positions in the representations, this
        function annotates whether the parents have a connector that connects
        them.
        
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
        
        Parameters
        ----------
        parent1_repres : list
            Representation of parent1 (aligned against parent2).
        parent2_repres : list
            Representation of parent2 (aligned against parent1).

        Returns
        -------
        connectors_table : 2D list
            This table stores at row i, column j the connector(s) available to
            link index i to index j.
        '''
        
        
        # XXX Code here to catch case when organisms were not set yet
        n = len(self.organism1)
        
        # 2D list where each item is an emtpy list
        connectors_table = [[ [] for i in range(n)] for j in range(n)]
        
        # Each parent representation is coupled with a tag ('p1' or 'p2')
        parents = [(self.organism1, 'p1'),
                   (self.organism2, 'p2')]
        
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









