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









