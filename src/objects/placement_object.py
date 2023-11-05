# -*- coding: utf-8 -*-
"""
Placement object

"""

import numpy as np

class PlacementObject:
    """
    Placement object
    
    """
    
    def __init__(self, organism_id, dna_sequence):
        """
        PlacementObject object constructor.

        Args:
            organism_id: str
            dna_sequence: str
        """
        
        self.organism_id = organism_id
        self.dna_sequence = dna_sequence
        
        # Initialize placement features
        
        self.energy = None
        self.recognizers_scores = []
        self.connectors_scores = []
        self.recognizers_positions = []
        self.connectors_positions = []
        self.recognizer_types = ""
    
    # Compile placement features
    
    def set_energy(self, energy):
        self.energy = energy
    
    def set_recognizer_scores(self, recog_scores):
        self.recognizers_scores = recog_scores
    
    def set_connectors_scores(self, conn_scores):
        self.connectors_scores = conn_scores
    
    def append_recognizer_position(self, recog_position):
        self.recognizers_positions.append(recog_position)
    
    def append_connector_position(self, connector_position):
        self.connectors_positions.append(connector_position)
    
    # Print placement (to standard output or to file)
    
    def print_placement_old(self, stdout=False, outfile=None):
        '''
        It prints out the placement as text.

        Parameters
        ----------
        stdout : bool, optional
            Whether to print to standard output the placement.
            The default is False.
        outfile : None or string, optional
            Whether to write to file the placement. If it's None, the placement
            is not written to file. Otherwise, it specifies the name of the
            outfile. The default is None.
        '''
        n = len(self.dna_sequence)
        recog_positions_line = ["-"] * n
        recog_scores_line = ["_"] * n
        conn_scores_line = ["_"] * n
        
        # Recognizers positions
        for i in range(len(self.recognizers_positions)):
            # start and stop DNA positions of recognizer i
            start, stop = self.recognizers_positions[i]
            
            # get the recognizer score and format it
            recog_score = self.recognizers_scores[i]
            recog_score_str = "{:.2f}".format(recog_score)
            
            # write recognizer placement
            for pos in range(start, stop):
                recog_positions_line[pos] = str(self.recognizer_types[i])
            
            # write recognizer score
            for c in range(len(recog_score_str)):
                if start + c < len(recog_scores_line):  # avoid going out of the seq
                    recog_scores_line[start + c] = recog_score_str[c]
        
        # Connectors positions
        for i in range(len(self.connectors_positions)):
            # start and stop DNA positions of connector i
            start, stop = self.connectors_positions[i]
            
            # get the connector score and format it
            connector_score = self.connectors_scores[i]
            connector_score_str = "{:.2f}".format(connector_score)
            
            gap_size = stop - start
            
            # if the gap is large, the connector score is written in the middle
            # if gap_size > len(connector_score_str) + 1:
            right_shift = int(np.ceil((gap_size - len(connector_score_str))/2))
            start += right_shift
            
            # # More centered printing for small gaps
            # else:
            #     start -= 2
            
            # write connector score
            for c in range(len(connector_score_str)):
                if start + c < len(conn_scores_line):  # avoid goin out of the seq
                    conn_scores_line[start + c] = connector_score_str[c]
        
        # print to stdout if required
        if stdout:
            print(self.dna_sequence)
            print("".join(recog_positions_line))
            print("".join(recog_scores_line))
            print("".join(conn_scores_line))
        
        # print to file if required
        if outfile != None:
            print(self.dna_sequence, file=outfile)
            print("".join(recog_positions_line), file=outfile)
            print("".join(recog_scores_line), file=outfile)
            print("".join(conn_scores_line), file=outfile)
    
    def print_placement(self, stdout=False, outfile=None):
        '''
        It prints out the placement as text.

        Parameters
        ----------
        stdout : bool, optional
            Whether to print to standard output the placement.
            The default is False.
        outfile : None or string, optional
            Whether to write to file the placement. If it's None, the placement
            is not written to file. Otherwise, it specifies the name of the
            outfile. The default is None.
        '''
        
        if not stdout and not outfile:
            return
        
        n = len(self.dna_sequence)
        recog_positions_line = [" "] * n
        recog_scores_line = [" "] * n
        conn_scores_line = ["_"] * n
        
        # Recognizers positions
        for i in range(len(self.recognizers_positions)):
            # start and stop DNA positions of recognizer i
            start, stop = self.recognizers_positions[i]
            
            # get the recognizer score and format it
            recog_score = self.recognizers_scores[i]
            recog_score_str = "{:.2f}".format(recog_score)
            if i == 0:
                conn_scores_line[:stop] = [" "] * stop
            if i == len(self.recognizers_positions) - 1:
                conn_scores_line[start:] = [" "] * (n - start)
            
            # write recognizer placement
            for pos in range(start, stop):
                recog_positions_line[pos] = str(self.recognizer_types[i])
                conn_scores_line[pos] = " "
            
            # write recognizer score
            for c in range(len(recog_score_str)):
                if start + c < len(recog_scores_line):  # avoid going out of the seq
                    recog_scores_line[start + c] = recog_score_str[c]
        
        # Connectors positions
        for i in range(len(self.connectors_positions)):
            # start and stop DNA positions of connector i
            start, stop = self.connectors_positions[i]
            
            # get the connector score and format it
            connector_score = self.connectors_scores[i]
            connector_score_str = "{:.2f}".format(connector_score)
            
            gap_size = stop - start
            
            # if the gap is large, the connector score is written in the middle
            # if gap_size > len(connector_score_str) + 1:
            right_shift = int(np.ceil((gap_size - len(connector_score_str))/2))
            start += right_shift
            
            # # More centered printing for small gaps
            # else:
            #     start -= 2
            
            # write connector score
            for c in range(len(connector_score_str)):
                if start + c < len(conn_scores_line):  # avoid goin out of the seq
                    conn_scores_line[start + c] = connector_score_str[c]
        
        # print to stdout if required
        if stdout:
            print("")
            print("".join(conn_scores_line))
            print("".join(recog_scores_line))
            print("".join(recog_positions_line))
            print(self.dna_sequence)
        
        # print to file if required
        if outfile != None:
            print("", file=outfile)
            print("".join(conn_scores_line), file=outfile)
            print("".join(recog_scores_line), file=outfile)
            print("".join(recog_positions_line), file=outfile)
            print(self.dna_sequence, file=outfile)
    
    def set_recognizer_types(self, rec_types):
        #print(rec_types)
        self.recognizer_types = rec_types




