# -*- coding: utf-8 -*-
"""
Organism object
It allocates the full data structure
"""

import numpy as np
from .placement_object import PlacementObject
from .shape_object import ShapeObject
from FMS.objects import shape_object as shape
from .connector_object import norm_cdf
from .connector_object import g
import _multiplacement
import decimal as dec
import json
from .connector_object import ConnectorObject
from .pssm_object import PssmObject
import scipy
import os
import pickle
from FMS.models import null

class ChainModel:
    """
    ChainModel object
    The ChainModel object essentially contains two vectors:
        - A vector of connector objects
        - A vector of recognizer objects
    The order of the elements in these vectors determine, implicitly, the
    connections between the elements.
    """

    _organism_counter = 0

    def __init__(self, _id: str) -> None:
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

        # self.is_precomputed = conf["PRECOMPUTE"]
        self.is_precomputed = True

        # !!! maximum length of PSSMs allowed
        # self.max_pssm_length = max_pssm_length

        # Map used by the placement algorithm
        # The list maps each row of the matrix of the placement scores onto a
        # column of a PSSM: each row is assigned a [pssm_idx, column_idx]
        self.row_to_pssm = []

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

    def to_json(self) -> list:
        """Get the chain model in JSON format
        """

        model = []
        for i in range(self.count_recognizers() - 1):
            model.append((self.recognizers[i]).to_json())
            model.append((self.connectors[i]).to_json())
        model.append((self.recognizers[-1]).to_json())
        return model
    
    def load_shape_null_models(self):
        """loads the serialzed null model dictionary into shape_object.null_models
        shape_object.null_models is a global variable that holds all null models for
        all shape recognizers for a given number of bins (specified by "NUM_BINS" 
        in the config.json file)

        Args:
            None

        Returns:
            None
        """
        if not os.path.isfile(".models/models"):
            shape.null_models = {}
            return
        
        if os.path.getsize(".models/models") == 0:
            shape.null_models = {}
            return
        
        with open(".models/models", "rb") as infile:
            try:
                shape.null_models = pickle.load(infile)
            except:
                shape.null_models = {}

            return

    def check_shape_null_models(self, n_bins):
        """Determines if null models must be computed or not, and calls 
        computing function if necessary.

        generates models for all shapes for lengths [5, MAX_LENGTH] if
        shape_object.null_models is empty, or empty for the specified number
        of intervals. If it isn't empty for the specified number of bins,
        and the maximum length that has been computed already is less than
        MAX_LENGTH, then it will compute the remaining null models and adds
        them to the dictionary. If all of the required null models have been
        computed already, it will do nothing and return None.

        Args:
            n_bins: the number of intervals in the null models for each shape

        Returns:
            None
        """

        if shape.null_models == {}:
            null.generate_range(5, 9 + 1,
                                shape.null_models, n_bins)
            return

        if n_bins not in shape.null_models.keys():
            null.generate_range(5, 9 + 1,
                                shape.null_models, n_bins)
            return

        if max(shape.null_models[n_bins]["mgw"].keys()) < 9:
            null.generate_range(max(shape.null_models[n_bins]["mgw"].keys()) + 1,
                                9 + 1,
                                shape.null_models, n_bins)
            return


    @staticmethod
    def import_model(file_name: str):
        """Import Chain Model from file

        Args:
            file_name: Name of the file with the chain model to read as an input

        Returns:
            a chain model objects read from the file
        """
        # global ShapeObject.null_models

        with open(file_name) as json_file:
            organism_json = json.load(json_file)[0]

        new_model = ChainModel(str(ChainModel._organism_counter))

        new_model.load_shape_null_models()
        num_bins = 25
        new_model.check_shape_null_models(num_bins)
        print(shape.null_models)
        shape.null_models = shape.null_models[num_bins]
        print(shape.null_models.keys())
        print(shape.null_models["helt"])

        new_org_recognizers = []  # The recognizers are collected here
        new_org_connectors = []  # The connectors are collected here

        for element in organism_json:
            if element["objectType"] == "pssm":
                new_org_recognizers.append(PssmObject(np.array(element["pwm"])))
            if element["objectType"] == "shape":
                new_org_recognizers.append(ShapeObject(element["recType"], element["length"], element["mu"], element["sigma"]))
            elif element["objectType"] == "connector":
                new_org_connectors.append(ConnectorObject(500, element["mu"], element["sigma"]))

        # Set recognizers and connectors of the organism
        new_model.set_recognizers(new_org_recognizers)
        new_model.set_connectors(new_org_connectors)
        new_model.set_row_to_pssm()

        # Check if the frequency values in the PWMs of the imported organisms
        # are consistent with  PSSM_NUM_OF_BINDING_SITES  set in the config file
        new_model.check_pwm_frequencies()
        new_model.flatten()

        ChainModel._organism_counter+=1

        return new_model
    
    def check_pwm_frequencies(self) -> None:
        '''
        Raises an error if the PWMs of the imported organisms are not compatible
        with the value of PSSM_NUM_OF_BINDING_SITES parameter specified in the
        config file.
        '''
        no_BSs = 100
        # no_BSs = config.configOrganismFactory["PSSM_NUM_OF_BINDING_SITES"]
        smallest_freq = dec.Decimal('1') / dec.Decimal(str(no_BSs))

        for rec_idx, rec in enumerate(self.recognizers):
            if rec.type != 'p':
                continue
            for pos in range(rec.length):
                for b in ['a','c','g','t']:
                    freq = rec.pwm[pos][b]
                    if dec.Decimal(str(freq)) % smallest_freq != 0:
                        raise Exception(
                            ("Imported organism has PWM frequencies that are not "
                                "compatible with the required number of binding sites "
                                "(PSSM_NUM_OF_BINDING_SITES parameter). The problem "
                                "occurred for the frequency of base " + b.upper() +
                                " at position " + str(pos) + " of the recognizer with "
                                "index " + str(rec_idx) + " of the organism "
                                "to be imported from the json file. "
                                "Indeed, " + str(freq) + " is not "
                                "k * " + str(smallest_freq) + " for any integer "
                                "value of k. "
                                "Please change the PSSM_NUM_OF_BINDING_SITES parameter "
                                "in the config file, or modify the organism to be "
                                "imported, accordingly to the desired number of "
                                "binding sites. ")
                        )

    def get_max_length(self, alpha):
        length = self.sum_recognizer_lengths
        print("length recognizers", length)
        mu = 0
        var = 0
        for connector in self.connectors:
            mu += connector._mu
            var += connector._sigma**2
        length += scipy.stats.norm(loc=mu, scale=var**0.5).ppf(alpha+(1-alpha)/2)
        return int(2*length)