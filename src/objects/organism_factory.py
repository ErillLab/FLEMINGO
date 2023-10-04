"""Organism Factory creates organisms, connectors and pssms
   It performs all operations dealing with multiple objects, 
   such as organism recombination.
"""

import random
import json
import numpy as np
from .shape_object import ShapeObject
from objects import shape_object
from models import null
from .organism_object import OrganismObject
from .connector_object import ConnectorObject
from .pssm_object import PssmObject
from .aligned_organisms_representation_object import AlignedOrganismsRepresentation
import copy
import decimal as dec
import pickle
import os


class OrganismFactory:
    """Factory class
    """

    def __init__(self, conf_org, conf_org_fac, conf_con, conf_pssm, p_rank,
                 conf_shape, min_seq_length, max_seq_length) -> None:
        """
        Instantiates an OrganismFactory object.
        Reads in the configuration paramaters for the factory and for all object
        types (organism, connector and PSSM recognizer). It also uses process
        rank (which is None in serial runs) to generate unique organism IDs when
        the program is run in parallel mode.
        """
        
        # Counter used to assign IDs to new organisms
        self._organism_counter = 0
        # Process rank (used to esure unique organism IDs when running in parallel mode)
        self._process_rank = p_rank
        # Length of shortest and longest DNA sequence
        self.min_seq_length = min_seq_length
        self.max_seq_length = max_seq_length
        # lambda parameter for Poisson distribution used to instantiate organisms.
        # lambda is the expected number of recognizers per organism
        self.avg_num_of_recognizers = conf_org_fac["AVG_NUM_OF_RECOGNIZERS"]
        
        self.probability_recombination = conf_org_fac["PROBABILITY_RECOMBINATION"]
        
        # MLE-based memetic drive
        if conf_org_fac["PERIODIC_MLE"] == None:
            self.periodic_mle = np.inf
        else:
            self.periodic_mle = conf_org_fac["PERIODIC_MLE"]
        self.probability_mle_pssm = conf_org_fac["PROBABILITY_MLE_PSSM"]
        self.probability_mle_shape = conf_org_fac["PROBABILITY_MLE_SHAPE"]
        self.probability_mle_connector = conf_org_fac["PROBABILITY_MLE_CONNECTOR"]
        self.probability_mle_pssm_to_shape = conf_org_fac["PROBABILITY_MLE_PSSM_TO_SHAPE"]
        self.probability_mle_shape_to_pssm = conf_org_fac["PROBABILITY_MLE_SHAPE_TO_PSSM"]
        
        # Minimum and maximum values allowed for connector's mu
        self.connector_min_mu = conf_org_fac["CONNECTOR_MIN_MU"]
        self.connector_max_mu = conf_org_fac["CONNECTOR_MAX_MU"]

        # Minimum and maximum values allowed for connector's sigma
        self.connector_min_sigma = conf_org_fac["CONNECTOR_MIN_SIGMA"]
        self.connector_max_sigma = conf_org_fac["CONNECTOR_MAX_SIGMA"]
        
        # Average length of PSSM
        self.recognizer_avg_length = conf_org_fac["RECOGNIZER_AVG_LENGTH"]
        
        # Number of binding sites used to generate the PWM
        self.pssm_number_of_binding_sites = conf_org_fac["PSSM_NUM_OF_BINDING_SITES"]
        self.pssm_vs_shape_probability = conf_org_fac["PSSM_VS_SHAPE_PROBABILITY"]
        
        # Assign organism, connector and pssm configurations
        self.conf_org = conf_org
        self.conf_con = conf_con
        self.conf_pssm = conf_pssm
        self.conf_shape = conf_shape
        self.load_shape_null_models()
        self.check_shape_null_models(conf_shape["NUM_BINS"])
        shape_object.null_models = shape_object.null_models[conf_shape["NUM_BINS"]]

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
        if not os.path.isfile("models/models"):
            shape_object.null_models = {}
            return
        
        if os.path.getsize("models/models") == 0:
            shape_object.null_models = {}
            return
        
        with open("models/models", "rb") as infile:
            try:
                shape_object.null_models = pickle.load(infile)
            except:
                shape_object.null_models = {}

            return

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

        if shape_object.null_models == {}:
            null.generate_range(5, self.conf_shape["MAX_LENGTH"] + 1,
                                shape_object.null_models, n_bins)
            return

        if n_bins not in shape_object.null_models.keys():
            null.generate_range(5, self.conf_shape["MAX_LENGTH"] + 1,
                                shape_object.null_models, n_bins)
            return

        if max(shape_object.null_models[n_bins]["mgw"].keys()) < self.conf_shape["MAX_LENGTH"]:
            null.generate_range(max(shape_object.null_models[n_bins]["mgw"].keys()) + 1,
                                self.conf_shape["MAX_LENGTH"] + 1,
                                shape_object.null_models, n_bins)
            return

        return

    def get_id(self) -> str:
        """ Returns a new unique organism ID (as a string). If the program is
        run in parallel the organism ID is composed of process rank and the
        organism counter ('r_n' format). For example, '3_19' would be the ID of
        the 19th organism generated by process 3. This ensures that the organism
        IDs are unique, without requiring MPI communication between processes.
        If the program is run in serial mode, the organism ID is simply the
        organism counter (as a string). For example, '21' would be the ID of the
        21th generated organism. """
        self._organism_counter += 1
        if self._process_rank is None:
            return str(self._organism_counter)
        else:
            return str(self._process_rank) + "_" + str(self._organism_counter)

    def get_organism(self) -> OrganismObject:
        '''
        It creates and returns a full organism datastructure.
        An organism contains essentially two lists:
            - a recognizer list
            - a connector list
        The index of these elements in the lists defines implicitly the
        connections between them.

        Returns:
            A new random organism based on JSON config file
        '''
        
        # instantiates organism with organism configuration and pssm columns
        new_organism = OrganismObject(self.get_id(), self.conf_org)
        
        '''
        The number of recognizer for the new organism is drawn from a Poisson
        distribution. The Poisson distribution is shifted so that the lowest
        bound is not 0, but `min_n_recogs`. The value of `min_n_recogs` is 1 unless
        a larger number is required to comply with the parameter that specifies
        the minimum number of nodes per organism in the config file. The shifted
        Poisson distribution is constructed so that the average number of nodes
        per organism matches the request from the config file.
        In other words, the Poisson distribution is used to decide how many nodes
        more than the minimum should an organism have, while ensureing the desired
        average number of nodes.
        However, if the number of recognizers drawn from the shifted Poisson
        exceeds the maximum allowed by the config, it is set to the maximum.
        '''
        # Draw from shifted Poisson
        n_recogs = self.from_shifted_poisson(self.avg_num_of_recognizers,
                                             new_organism.min_num_of_recognizers,
                                             new_organism.max_num_of_recognizers)
        
        # for each recognizer in the organism except for the last
        for i in range(n_recogs - 1):
            # Add one recognizer
            new_recognizer = self.create_recognizer(self.recognizer_avg_length)
            new_organism.recognizers.append(new_recognizer)
            # Add one connector
            # _mu = random.randint(self.connector_min_mu, self.connector_max_mu)
            # _sigma = random.randint(self.connector_min_sigma, self.max_sigma)
            # new_connector = ConnectorObject(self.conf_con, self.max_seq_length, _mu, _sigma)
            new_connector = self.create_connector()
            new_organism.connectors.append(new_connector)
        # Add last recognizer
        new_recognizer = self.create_recognizer(self.recognizer_avg_length)
        new_organism.recognizers.append(new_recognizer)
        # Set attribute that will map organism nodes to alignment matrix rows
        new_organism.set_row_to_pssm()
        new_organism.flatten()
        return new_organism
    
    def from_shifted_poisson(self, average, minimum, maximum):
        '''
        Randomly draws an integer from a Poisson distribution that is shifted
        so that the smallest number returned is the number specified by the
        parameter `minimum`. It also ensures that the expected value is the
        number specified by the parameter `average`. If the randomly drawn
        number is larger than `maximum` it is set to `maximum`.
        '''
        # Draw from shifted Poisson
        n = np.random.poisson(average - minimum)
        n += minimum
        # Check upper bound (if specified)
        if maximum != None:
            n = min(n, maximum)
        return n
    
    def create_connector(self) -> ConnectorObject:
        '''
        It returns a connector object with its internal parameters (mu, sigma)
        assigned.
        '''
        return ConnectorObject(self.conf_con, self.max_seq_length)

    def create_recognizer(self, length=None):
        '''
        Randomly creates either a PSSM or shape recognizer.

        Args:
            length:
                length of the recognizer to be created. Default is None.
                Passed to either create_pssm or create_shape.

        Returns:
            A recognizer with the specified length (if its valid), or a
            randommly chosen length if not specified.
        '''
        if random.random() < self.pssm_vs_shape_probability:
            return self.create_pssm(length)
        else:
            return self.create_shape(length)

    def create_pssm(self, length=None) -> PssmObject:
        '''
        It return a PSSM object. The length can be specified.

        Args:
            length:
                length (number of columns) of the PSSM. Default is None.
                If None, the length is randomly drawn from a modified Poisson
                distribution that respects the required average and bounds.
        
        Returns:
            A pssm object with an initializated PWM
        '''
        if length == None:
            # length = self.recognizer_avg_length
            length = self.from_shifted_poisson(self.recognizer_avg_length,
                                               self.conf_pssm["MIN_LENGTH"],
                                               self.conf_pssm["MAX_LENGTH"])
        
        # Generate as many PSSM columns as needed
        pwm = []
        for _ in range(length):
            pwm.append(self.get_pwm_column())

        return PssmObject(np.array(pwm), self.conf_pssm)
    
    def create_shape(self, length=None, rec_type=None) -> ShapeObject:
        '''
        It return a shape object. The length can be specified (>= 5). The shape
        type can be specified (one of {"mgw", "prot", "roll", "helt"}).

        Args:
            length:
                length of the shape recognizer. Default is None.
                If None or less than 5, then 5 is used.
            rec_type:
                type of shape recognizer. Default is None.
                If None, a random type is chosen.
        
        Returns:
            A shape object
        '''
        if length == None or length < self.conf_shape["MIN_LENGTH"]:
            length = self.conf_shape["MIN_LENGTH"]
        if length > self.conf_shape["MAX_LENGTH"]:
            length = self.conf_shape["MAX_LENGTH"]
        
        if rec_type == None:
            rec_type = random.choice(['mgw', 'prot', 'roll', 'helt'])
        return ShapeObject(rec_type, length, self.conf_shape)
    
    def get_pwm_column(self) -> dict:
        """Generates a single column for a PWM

        Returns:
            a random probability for each base [a, c, g, t]
        """
        
        number_of_BS_to_assign = self.pssm_number_of_binding_sites
        
        # we assign a number of sites (a count) to each row (to each base)
        counts: list = []
        
        # First three count values
        for item in range(3,0,-1):
            # make sure that no row gets less than one site
            if number_of_BS_to_assign > item:
                count = random.randint(1, number_of_BS_to_assign-item)
                number_of_BS_to_assign -= count
            else:
                count = 1
                number_of_BS_to_assign -= count
            counts.append(count)
        
        # All the remaining unassigned counts go to the last base
        counts.append(number_of_BS_to_assign)

        # Shuffle list so that highest count isn't always on first base ("A")
        random.shuffle(counts)

        # Convert counts to probabilities
        np_probabilities = np.array(counts) / self.pssm_number_of_binding_sites
        probabilities = np_probabilities.tolist()
        
        return {
            "a": probabilities[0],
            "g": probabilities[1],
            "c": probabilities[2],
            "t": probabilities[3],
        }
    
    def import_organisms(self, file_name: str) -> list:
        """Import Organisms from file

        Args:
            file_name: Name of the file with the organisms to read as an input

        Returns:
            a list of organisms objects read from the file
        """
        organism_list = []

        with open(file_name) as json_file:
            organism_json = json.load(json_file)

        for organism in organism_json:

            new_organism = OrganismObject(
                self.get_id(), self.conf_org
            )
            
            new_org_recognizers = []  # The recognizers are collected here
            new_org_connectors = []  # The connectors are collected here
            
            for element in organism:
                if element["objectType"] == "pssm":
                    new_org_recognizers.append(self.import_pssm(element))
                if element["objectType"] == "shape":
                    new_org_recognizers.append(self.import_shape(element))
                elif element["objectType"] == "connector":
                    new_org_connectors.append(self.import_connector(element))
            
            # Set recognizers and connectors of the organism
            new_organism.set_recognizers(new_org_recognizers)
            new_organism.set_connectors(new_org_connectors)
            new_organism.set_row_to_pssm()
            organism_list.append(new_organism)
        
        # Check if the frequency values in the PWMs of the imported organisms
        # are consistent with  PSSM_NUM_OF_BINDING_SITES  set in the config file
        self.check_pwm_frequencies_of_imported_organisms(organism_list)
        for organism in organism_list:
            organism.flatten()
        return organism_list
    
    def check_pwm_frequencies_of_imported_organisms(self, imported_organisms: list):
        '''
        Raises an error if the PWMs of the imported organisms are not compatible
        with the value of PSSM_NUM_OF_BINDING_SITES parameter specified in the
        config file.
        '''
        no_BSs = self.pssm_number_of_binding_sites
        smallest_freq = dec.Decimal('1') / dec.Decimal(str(no_BSs))
        
        for org_idx, org in enumerate(imported_organisms):
            for rec_idx, rec in enumerate(org.recognizers):
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
                                 "index " + str(rec_idx) + " of the organism with "
                                 "index " + str(org_idx) + " in the list of organisms "
                                 "to be imported from the json file. "
                                 "Indeed, " + str(freq) + " is not "
                                 "k * " + str(smallest_freq) + " for any integer "
                                 "value of k. "
                                 "Please change the PSSM_NUM_OF_BINDING_SITES parameter "
                                 "in the config file, or modify the organism to be "
                                 "imported, accordingly to the desired number of "
                                 "binding sites. ")
                            )
    
    def import_connector(self, connector: dict) -> ConnectorObject:
        """Import Connector from JSON object

        Args:
            connector: connector in dictionary format

        Returns:
            Connector object from given connector dictionary
        """
        return ConnectorObject(self.conf_con, self.max_seq_length, connector["mu"], connector["sigma"])


    def import_pssm(self, pssm: dict) -> PssmObject:
        """Import PSSM from JSON object

        Args:
            pssm: pssm recognizer in dictionary format

        Returns:
            PSSM Object from given  pssm dictionary
        """
        return PssmObject(np.array(pssm["pwm"]), self.conf_pssm)

    def import_shape(self, shape: dict) -> ShapeObject:
        """Import shape recognizer from JSON object

        Args:
            shape: shape recognizer in dictionary format

        Returns:
            Shape Object from given shape dictionary
        """
        return ShapeObject(shape["recType"], shape["length"], self.conf_shape, shape["mu"], shape["sigma"])

    def export_organisms(self, a_organisms: list, filename: str) -> None:
        """Export a list of organisms to JSON format

        Args:
            a_organisms: list of organisms to export
            filename: name of the file to export all the organisms
        """
        
        list_json_organisms = []
        for o_organism in a_organisms:
            organism = []
            for i in range(o_organism.count_recognizers() - 1):
                organism.append(self.export_recognizer(o_organism.recognizers[i]))
                organism.append(self.export_connector(o_organism.connectors[i]))
            organism.append(self.export_recognizer(o_organism.recognizers[-1]))
            list_json_organisms.append(organism)
        
        with open(filename, "w+") as json_file:
            json.dump(list_json_organisms, json_file, indent=2)

    def export_connector(self, o_connector: ConnectorObject) -> dict:
        """Export connector object

        Args:
            o_connector: Connector object to export

        Returns:
            Connector in dictionary format
        """
        connector = {}
        connector["objectType"] = "connector"
        connector["mu"] = o_connector._mu
        connector["sigma"] = o_connector._sigma

        return connector

    def export_recognizer(self, o_rec) -> dict:
        """calls the appropriate export function depending on
        the type of the given recognizer

        Args:
            o_rec: Recognizer object to export

        Returns:
            Recognizer in dictionary format
        """
        if o_rec.get_type() == 'p':
            return self.export_pssm(o_rec)
        else:
            return self.export_shape(o_rec)
        
    def export_pssm(self, o_pssm: PssmObject) -> dict:
        """Export PSSM object

        Args:
            o_pssm: PSSM object to export

        Returns:
            pssm in dictionary format
        """
        pssm = {}
        pssm["objectType"] = "pssm"
        pssm["pwm"] = o_pssm.pwm.tolist()
        return pssm
    
    
    def export_shape(self, o_shape: ShapeObject) -> dict:
        """Export Shape object

        Args:
            o_shape: Shape object to export

        Returns:
            shape in dictionary format
        """
        shape = {}
        shape["objectType"] = "shape"
        shape["recType"] = o_shape.type
        shape["mu"] = o_shape._mu
        shape["sigma"] = o_shape._sigma
        shape["length"] = o_shape.length
        return shape

    def get_children(self, par1, par2, reference_dna_seq, pos_dna_sample):
        '''
        Implements the recombination operator.
        Fisrt, an abstract representation of the aligned parents is produced.
        Then, the symbols of these representations are swapped, generating the
        representations of the children. In this way we define which nodes will
        go in each child.
        Finally, the actual children organism objects are compiled accordingly.
        
        Parameters
        ----------
        par1 : OrganismObject
            First parent.
        par2 : OrganismObject
            Second parent.
        reference_dna_seq : string
            The DNA sequence on which the parents are placed, so that they
            become 'aligned' one against the other.
        pos_dna_sample : list
            A random subset of the positive set. It's used in case it's
            necessary to make some synthetic connectors (for the children). The
            average distance and the standard deviation will be estimated using
            this sample.

        Returns
        -------
        child1 : OrganismObject
            First child.
        child2 : OrganismObject
            Second child.
        '''
        
        # Initialize child 1 as an empty organism
        child1 = OrganismObject(self.get_id(), self.conf_org)
        
        # Initialize child 2 as an empty organism
        child2 = OrganismObject(self.get_id(), self.conf_org)
        
        # Place the parents on all the sequences in the sample of the positive set
        par1_placements, par2_placements = self.store_parents_placemnts(par1, par2, pos_dna_sample)
        
        # Representation of the two parents aligned
        parents_repres = AlignedOrganismsRepresentation(par1, par2, reference_dna_seq)
        # Table storing info about what connectors are available to cover the possible spans
        parents_repres.annotate_available_connectors()
        
        # Representation of the two recombined children aligned
        children_repres = self.get_aligned_children_repr(parents_repres, child1._id, child2._id)
        
        # Assemble child 1
        # Write the assembly instructions
        child1.set_assembly_instructions(children_repres, 'org1')
        
        # Compile child 1
        self.compile_recognizers(child1, par1, par2)
        self.compile_connectors(child1, par1, par2, parents_repres,
                                par1_placements, par2_placements)
        
        # Assemble child 2
        # Write the assembly instructions
        child2.set_assembly_instructions(children_repres, 'org2')
        
        # Compile child 2
        self.compile_recognizers(child2, par1, par2)
        self.compile_connectors(child2, par1, par2, parents_repres,
                                par1_placements, par2_placements)
        child1.flatten()
        child2.flatten()
        return child1, child2
    
    def store_parents_placemnts(self, parent1, parent2, dna_seq_set):
        '''
        Places each parent on each DNA sequence in the given list of sequences.
        Returns all the placements in a list, for each organism.
        '''
        p1_placements = []
        p2_placements = []
        
        for dna_seq in dna_seq_set:
            p1_placements.append(parent1.get_placement(dna_seq))
            p2_placements.append(parent2.get_placement(dna_seq))
        
        return p1_placements, p2_placements
    
    def get_aligned_children_repr(self, parents_repres, child1_id, child2_id):
        '''
        This function swaps parts of the representations of the two parents, in
        order to get the representations of the two children.
        '''
        
        # Initialize representation of the children as identical copies of the parents
        children_repres = copy.deepcopy(parents_repres)
        children_repres.set_children_IDs(child1_id, child2_id)
        
        # Within each unit, perform a swap with 50% probability
        for (start, stop) in children_repres.units:
            if random.random() < 0.5:
                # Perform the swapping, which means that the part from parent1 will
                # end up into child2, and the part from parent2 will end up into child1
                children_repres.swap_unit(start, stop)
        
        return children_repres
    
    def compile_recognizers(self, child_obj, parent1, parent2):
        '''
        It appends to the given organism (child_obj) the required recognizers
        from the two parents, in the right order, according to the
        assembly_instructions attribute of the organism.
        '''
        
        for recog_name in child_obj.assembly_instructions['recognizers']:
            parent, recog_idx = recog_name.split('_')
            
            if parent == 'p1':
                recog = copy.deepcopy(parent1.recognizers[int(recog_idx)])
            elif parent == 'p2':
                recog = copy.deepcopy(parent2.recognizers[int(recog_idx)])
            
            # Add recognizer to organism
            child_obj.append_recognizer(recog)
    
    def compile_connectors(self, child_obj, par1, par2, parents_repr,
                           par1_placements, par2_placements):
        '''
        It appends to the given organism (child_obj) the required connectors
        from the two parents, in the right order, according to the
        assembly_instructions attribute of the organism.
        
        When a required connector was not available in the parents,
        assembly_instructions requires to synthesize a new connector. This is
        done by calling the make_synthetic_connector method.
        '''
        # Read assembly instructions
        recogs_names = child_obj.assembly_instructions['recognizers']
        connectors_names = child_obj.assembly_instructions['connectors']
        conn_adj_vals = child_obj.assembly_instructions['connectors_adjustments']
        
        for i in range(len(connectors_names)):
            connector_name = connectors_names[i]
            
            # "synth" means that the connector needs to be synthesized, because
            # it was not available in any of the two parents
            if connector_name[:5] == 'synth':                
                left_idx, right_idx = connector_name.split('_')[1:]
                
                # mu and sigma will be estimated for the gap between a left and a
                # right recognizers (on a small sample from the positive dataset).
                # Chose the left and right recognizers
                p1_left = parents_repr.organism1[int(left_idx)]
                p2_left = parents_repr.organism2[int(left_idx)]
                p1_right = parents_repr.organism1[int(right_idx)]
                p2_right = parents_repr.organism2[int(right_idx)]
                
                # When possible, chose them from the same parent
                if p1_left != '-' and p1_right != '-':
                    recog_L_name, recog_R_name = p1_left, p1_right
                elif p2_left != '-' and p2_right != '-':
                    recog_L_name, recog_R_name = p2_left, p2_right
                else:
                    if p1_left != '-':
                        recog_L_name, recog_R_name = p1_left, p2_right
                    else:
                        recog_L_name, recog_R_name = p2_left, p1_right
                # Make an appropriate connector
                conn = self.make_synthetic_connector(
                    recog_L_name, recog_R_name, par1_placements, par2_placements)
            
            # Else, the connector can be grabbed from one of the parents
            else:
                parent, conn_idx = connector_name.split('_')
                conn_idx = int(conn_idx)
                
                # Copy connector
                if parent == 'p1':
                    # Re-use connector from parent 1
                    conn = copy.deepcopy(par1.connectors[conn_idx])
                elif parent == 'p2':
                    # Re-use connector from parent 2
                    conn = copy.deepcopy(par2.connectors[conn_idx])
                
                # Apply left/right adjustments to the value of mu
                adjusted = False
                
                # Left and Right recognizers
                l_rec_name = recogs_names[i]
                r_rec_name = recogs_names[i + 1]
                
                # Apply LEFT adjustment if necessary
                if l_rec_name.split("_")[0] != parent:
                    left_adj, right_adj = conn_adj_vals[parent][conn_idx]
                    if left_adj != 0:
                        conn._mu += left_adj
                        conn.apply_mu_bounds()  # Avoide negative mu or too large mu
                        adjusted = True
                # Apply RIGHT adjustment if necessary
                if r_rec_name.split("_")[0] != parent:
                    left_adj, right_adj = conn_adj_vals[parent][conn_idx]
                    if right_adj != 0:
                        conn._mu += right_adj
                        conn.apply_mu_bounds()  # Avoide negative mu or too large mu
                        adjusted = True
                # If mu was adjusted, recompute the scores of PDF and CDF values
                if adjusted:
                    conn.set_precomputed_pdfs_cdfs()
            
            # Add connector to organism
            child_obj.append_connector(conn)
    
    def make_synthetic_connector(self, recog_left_name, recog_right_name,
                                 p1_placements, p2_placements):
        '''
        This function is used to generate an appropriate connector to link two
        recognizers of a child, when no one of the connectors from the parents
        is appropriate.

        Parameters
        ----------
        recog_left_name : string
            It identifies the recognizer to the left.
        recog_right_name : string
            It identifies the recognizer to the right.
        p1_placements : list
            List of placements for parent 1.
            Note:
                Placements are assumed to be paired. I.e., p1_placements[i] and
                p2_placements[i] are placements on the same DNA sequence
        p2_placements : list
            List of placements for parent 2.
            Note:
                Placements are assumed to be paired. I.e., p1_placements[i] and
                p2_placements[i] are placements on the same DNA sequence
        
        Returns
        -------
        synthetic_connector : ConnectorObject
            A new connector that connects the two input recognizers from two
            different parents, for which we didn't have an available connector.
            mu and sigma are inferred from the input placements. There are n
            placements of parent 1 on n DNA sequences and n placements of
            parent 2 on those same n DNA sequences. Therefore, for each of the
            n DNA sequences, we can estimate the relative distance of those two
            connectors. These stats are used to guess a reasonable value for mu
            and sigma.
        '''
        # Get parent and node index that specify the left recognizer
        recog_left_parent, recog_left_idx = recog_left_name.split('_')
        recog_left_idx = int(recog_left_idx)
        # Get parent and node index that specify the right recognizer
        recog_right_parent, recog_right_idx = recog_right_name.split('_')
        recog_right_idx = int(recog_right_idx)
        
        n_DNA_sequences = len(p1_placements)
        gap_values = []
        for i in range(n_DNA_sequences):
            
            # Placement of the left recognizer
            if recog_left_parent == 'p1':
                placement = p1_placements[i]
            elif recog_left_parent == 'p2':
                placement = p2_placements[i]
            L_first_bp, L_last_bp = self.get_recog_pos_on_DNA_seq(placement, recog_left_idx)
            
            # Placement of the right recognizer
            if recog_right_parent == 'p1':
                placement = p1_placements[i]
            elif recog_right_parent == 'p2':
                placement = p2_placements[i]
            R_first_bp, R_last_bp = self.get_recog_pos_on_DNA_seq(placement, recog_right_idx)
            
            # Gap between the left and the right recognizers
            distance = R_first_bp - L_last_bp
            gap = distance - 1
            gap_values.append(gap)
        
        # ESTIMATE SIGMA
        # Compute sample standard deviation as an estimate for sigma
        sigma = np.std(gap_values, ddof=1)
        # Avoid setting sigma to 0
        if sigma < 0.1:  # !!! temporarily hard-coded lower-bound
            sigma = 0.1
        
        # ESTIMATE MU
        '''
        mu could be estimated as the average gap. However, due to the typically
        small sample size, and the fact that in most use cases we expect several
        gaps to be optimal and others to be random values (coming from bad
        placements), we set mu to be equal to one of the gap values, instead of
        the average. This will sometimes lead to bad mu estimates, but often
        enough it will guess a clean estimate mu, avoiding the effect of bad
        placements. Another approach to this problem could have been to use the
        median instead of the average, but "bad placements" can produce
        consistently small/consistently large gaps, so the median could easily
        fail. Therefore we just choose one gap value to be our guess for mu.
        '''
        mu = gap_values[0]
        # Avoid negative mu
        mu = max(0, mu)
        
        # Return the synthetic connector
        return ConnectorObject(self.conf_con, self.max_seq_length, mu, sigma)
    
    def get_recog_pos_on_DNA_seq(self, org_placement, recog_idx):
        '''
        For a given placement and a given node index, it returns the first and
        last bp positions occupied by that node.

        Parameters
        ----------
        org_placement : dictionary
            The placement of an organism on a DNA sequence.
        node_idx : int
            Index of the desired node.

        Returns
        -------
        first_bp : int
            First DNA position occupied by the specified node.
        last_bp : int
            Last DNA position occupied by the specified node.
        '''
        
        start, stop = org_placement.recognizers_positions[recog_idx]
        
        first_bp = start  # First bp occupied
        last_bp = stop -1  # Last bp occupied
        
        return (first_bp, last_bp)
    
    # =========================================================================
    # XXX  Code for MLE
    # =========================================================================
    
    def clone_organism(self, organism):
        '''
        Returns an organism with the same parameters as the input organism (but
        different ID). The node tags in the 'assembly_instructions' attribute
        report that every node comes from parent 1 ('p1_') because there is
        only one parent (which is being cloned).
        '''
        clone = copy.deepcopy(organism)
        # Assign ID
        clone.set_id(self.get_id())
        # Clone assembly instructions
        clone.assembly_instructions['p1'] = organism._id
        clone.assembly_instructions['p2'] = None
        clone.assembly_instructions['recognizers'] = "cloned from " + organism._id
        clone.assembly_instructions['connectors'] = "cloned from " + organism._id
        return clone
    
    def mle_connector(self, conn_idx, placements):
        '''
        Given an organism that shows the placements specified by the parameter
        `placements`, it analyses the placements of the connector specified by
        the index `conn_idx`.
        It returns a new connector whose parameters are obtained through maximum
        likelihood estimation (MLE), to optimize binding energy on the positive
        set, given the observed placements on the positive set of the original
        connector.
    
        Parameters
        ----------
        conn_idx : int
            Index of the connector to be optimized according to MLE.
        placements : PlacementObject
            List of placements of the organism being optimized according to MLE.
    
        Returns
        -------
        conn : ConnectorObject
            A new connector object. Its parameters optimize binding energy on the
            positive set and were found based on MLE.
        '''
        
        # Observed gaps in the positive set
        gaps = []
        for plc in placements:
            start, stop = plc.connectors_positions[conn_idx]
            gaps.append(stop - start)
        
        # gap stats
        mu = np.mean(gaps)
        sigma = np.std(gaps, ddof=1)  # !!! ddof: 1 or 0? (sample VS population)
        
        # Create energy-optimized connector
        conn = self.create_connector()
        conn.set_mu(mu)
        conn.set_sigma(sigma)
        conn.set_precomputed_pdfs_cdfs()
        return conn
    
    def mle_pssm(self, recog_idx, placements):
        '''
        Given an organism that shows the placements specified by the parameter
        `placements`, it analyses the placements of the PSSM-recognizer
        specified by the index `recog_idx`.
        It returns a new PSSM whose parameters are obtained through maximum
        likelihood estimation (MLE), to optimize binding energy on the positive
        set, given the observed placements on the positive set of the original
        PSSM.
    
        Parameters
        ----------
        recog_idx : int
            Index of the PSSM to be optimized according to MLE.
        placements : PlacementObject
            List of placements of the organism being optimized according to MLE.
        organism : OrganismObject
            Organism the recognizer belongs to.
    
        Returns
        -------
        PssmObject :
            A new PSSM. Its parameters optimize binding energy on the
            positive set and were found based on MLE.
        '''
        
        # Recognized k-mers from the positive set
        instances = []
        for plc in placements:
            start, stop = plc.recognizers_positions[recog_idx]
            instances.append(plc.dna_sequence[start:stop].lower())
        
        # Instances matrix
        list_2D = []
        for inst in instances:
            list_2D.append([nucl for nucl in inst])
        inst_mat = np.array(list_2D)
        
        # PWM
        pwm = []
        pssm_len = len(instances[0])
        for i in range(pssm_len):
            column = inst_mat[:,i]
            # Nucleotides in `unique` and corresponding counts in `counts`
            unique, counts = np.unique(column, return_counts=True)
            # Report absent nucleotides as elements with 0 counts
            for nucl in ['a','c','g','t']:
                if nucl not in unique:
                    unique = np.append(unique, nucl)
                    counts = np.append(counts,0)
            # Observed frequencies
            frequencies = counts / sum(counts)
            # Convert observed frequencies into frequencies that are compatible
            # with the assumed number of binding sites (specified in the config)
            compatible_freqs = self.get_compatible_freqs(frequencies)
            # Compile PSSM column with compatible frequencies
            column_dict = dict(zip(unique, compatible_freqs))
            pwm.append(column_dict)
        
        # return PSSM object
        return PssmObject(np.array(pwm), self.conf_pssm)
    
    def get_compatible_freqs(self, column_freqs):
        '''
        This function turns the four given frequencies (in `column_freqs`) into
        frequencies that can be used in PWMs models of objects of the PssmObject
        class. Indeed, the frequencies of a PssmObject are constrained to be
        compatible with an assumed number of binding sites (the total number of
        counts), N.
        
        For example, if N=100, 0.01 or 0.53 are acceptable frequencies, while
        0.000013 or 0.4278 are not.
        
        First it assigns compatible counts with the 'get_compatible_counts'
        function. Then it turns those counts back to frequencies. When doing so,
        it makes sure that no approximation or numerical error turns those
        frequencies into numbers with extra decimal values (this is achieved
        using the 'decimal' library). The method ensures that frequencies are
        stored accurately (with the correct number of decimal places) and that
        their sum is exactly 1.
        '''
        # Compatible counts
        compatible_counts = self.get_compatible_counts(column_freqs)
        # Number of instances
        num_inst = sum(compatible_counts)
        if num_inst != self.pssm_number_of_binding_sites:
            raise ValueError(
                "Something went wrong. Number of instances is {}, but should have been {}.".format(
                    num_inst, self.pssm_number_of_binding_sites))
        # Back to frequencies
        compatible_freqs = compatible_counts / num_inst
        
        # To avoid extra decimals due to floating point errors:
        # Define the number of decimal digits actually required.
        # All frequencies are multiples of 1/N, where N is the required number of
        # binding sites (from the config file)
        smallest_freq = dec.Decimal('1') / dec.Decimal(str(num_inst))
        # Number of decimals in the smallest frequency
        num_decimals = len(str(smallest_freq).split(".")[1])
        # No frequency value needs more decimal digits than  smallest_freq.
        # Therefore we can round according to  num_decimals
        compatible_freqs = np.round(compatible_freqs, num_decimals)
        
        if float(sum([dec.Decimal(str(f)) for f in compatible_freqs])) != 1.0:
            raise ValueError(("Sum of freqs is not 1. Freqs are {}".format(compatible_freqs)))
        
        return compatible_freqs
    
    def get_compatible_counts(self, column_freqs):
        '''
        It assigns a count (an integer number) to each of the four input
        frequencies (the `column_freqs` argument) as faithfully as possible,
        while ensuring that the sum of the counts is N, where N is defined by
        the PSSM_NUM_OF_BINDING_SITES parameter in the config file.
        
        The four counts are assigned according to the "Largest Remainder Method"
        (LRM). This algorithm rounds down the scaled frequencies to the closest
        smaller integer. The remainder is distributed prioritizing the values
        that were rounded down the most. This is a commonly applied method to
        turn frequencies into counts while minimizing errors (minimizing the
        distance between the frequencies based on the output counts and the
        original frequencies).
        
        However, only N-4 counts are assigned proportionally to the input
        frequencies (according to the LRM). One count per nucleotide is
        assigned by defalut, serving as a pseudocount.
        
        RATIONALE:
            This function is useful to turn given frequencies into frequencies
            that can be used in the PWMs models by objects of the PssmObject
            class. Indeed, the frequencies of a PssmObject are constrained to
            be compatible with an assumed number of binding sites (the total
            number of counts), N. For example, if N=100, 0.01 or 0.53 are
            acceptable frequencies, while 0.000013 or 0.4278 are not.
        '''
        
        # Apply Largest Remainder Method to get a number of counts for each base
        
        # Four counts are not assigned yet because one count per base will be
        # added as a pseudocount
        n_counts_to_assign = self.pssm_number_of_binding_sites - 4
        non_int_counts = column_freqs * n_counts_to_assign
        int_counts = np.floor(non_int_counts).astype(int)
        order = np.flip(np.argsort(non_int_counts - int_counts))
        remainder = n_counts_to_assign - sum(int_counts)
        for i in range(remainder):
            int_counts[order[i]] += 1
        
        # Now we can add one count per base (the pseudocounts) so that the total
        # number of assigned counts will be equal to `self.pssm_number_of_binding_sites`
        int_counts = int_counts + 1
        return int_counts
    
    def get_shape_score(self, sequence, shape_type):
        '''
        Returns the shape value of the k-mer specified by the parameter `sequence`.
        The returned shape value is the one for the type of shape specified by
        the parameter `shape_type`.
        '''
        if shape_type == 'mgw':
            return null.get_sequence_mgw_score(sequence)
        elif shape_type == 'prot':
            return null.get_sequence_prot_score(sequence)
        elif shape_type == 'roll':
            return null.get_sequence_roll_score(sequence)
        elif shape_type == 'helt':
            return null.get_sequence_helt_score(sequence)
        else:
            raise ValueError('Unknown shape_type. shape_type should be "mgw" / ' +
                             '"prot" / "roll" / "helt".')
    
    def mle_shape(self, recog_idx, placements, shape_type):
        '''
        Given an organism that shows the placements specified by the parameter
        `placements`, it analyses the placements of the PSSM-recognizer
        specified by the index `recog_idx`.
        It returns a new shape-recognizer whose parameters are obtained through
        maximum likelihood estimation (MLE), to optimize binding energy on the
        positive set, given the observed placements on the positive set of the
        original shape-recognizer.
    
        Parameters
        ----------
        recog_idx : int
            Index of the shape-recognizer to be optimized according to MLE.
        placements : PlacementObject
            List of placements of the organism being optimized according to MLE.
        shape_type : str
            Type of DNA shape recognizer ("mgw"/"helt"/"roll"/"prot").
    
        Returns
        -------
        ShapeObject :
            A new shape-recognizer. Its parameters optimize binding energy on
            the positive set and were found based on MLE.
        '''
        # Recognized k-mers from the positive set
        instances = []
        for plc in placements:
            start, stop = plc.recognizers_positions[recog_idx]
            instances.append(plc.dna_sequence[start:stop].lower())
        
        # Shape length
        shape_length = len(instances[0])
        
        # observed shape values
        scores = []
        for sequence in instances:
            scores.append(self.get_shape_score(sequence, shape_type))
        
        # shape value stats
        mu = np.mean(scores)
        sigma = np.std(scores, ddof=1)  # !!! ddof: 1 or 0? (sample VS population)
        
        # Create energy-optimized shape-recognizer
        shape = self.create_shape(shape_length, shape_type)
        shape._mu = mu
        shape._sigma = sigma
        shape.set_alt_model()
        return shape
    
    def mle_org(self, organism, placements):
        '''
        Given an organism (`organism` argument) and its placements (`placements`
        argument), it returns a new organism, where the parameters of the nodes
        are obtained through maximum likelihood estimation (MLE). Each type of
        node has a certain probability of being optimized through MLE, which
        means that its binding instances are collected and the parameters are
        firectly inferred from that collection.
        For recognizers, there's also the possibility of using the observed
        instances of a PSSM to obtain an MLE-optimized Shape recognizer (of any
        type), or using the observed instances of a Shape recognizer to obtain
        an MLE-optimized PSSM.
        
        The probabilities of such events are the following:
            
            probability_mle_connector :
                probability that a connector undergoes MLE.
            probability_mle_pssm :
                probability that a PSSM undergoes MLE.
            probability_mle_shape :
                probability that a Shape recognizer undergoes MLE.
            probability_mle_pssm_to_shape :
                probability that a PSSM is turned into a MLE-optimized Shape
                recognizer. The shape type is chosen randomly.
            probability_mle_shape_to_pssm :
                probability that a Shape recognizer is turned into a
                MLE-optimized PSSM.
        
        Parameters
        ----------
        organism : OrganismObject
            Organism whose parameters need to be optimized according to MLE.
        placements : list
            List of placements (on the positive set) of the organism to be
            optimized according to MLE. Each placement is a PlacementObject
        
        Returns
        -------
        OrganismObject :
            A new organism that has some nodes optimized through MLE.
        '''
        # Optimize connectors
        mle_connectors = []
        for i in range(organism.count_connectors()):
            # MLE-optimized connector
            # [with probability `probability_mle_connector`]
            if random.random() < self.probability_mle_connector:
                mle_connectors.append(self.mle_connector(i, placements))
            # Leave the connector unchanged
            # [with probability 1 - probability_mle_connector]
            else:
                mle_connectors.append(copy.deepcopy(organism.connectors[i]))
        
        # Optimize recognizers
        mle_recognizers = []
        for i in range(organism.count_recognizers()):
            
            # PSSM
            if organism.recognizers[i].is_pssm():
                x = random.random()
                
                # (I) MLE-optimized PSSM
                # [with probability `probability_mle_pssm`]
                if x < self.probability_mle_pssm:
                    mle_recognizers.append(self.mle_pssm(i, placements))
                
                # (II) Turn the PSSM into a Shape (MLE-optimized)
                # [with probability `probability_mle_pssm_to_shape`]
                elif (x - self.probability_mle_pssm) < self.probability_mle_pssm_to_shape:
                    # Check that the PSSM is long enough to become a shape
                    if organism.recognizers[i].length >= 5:
                        # The type of shape is chosen randomly
                        shape_type = random.choice(['mgw','prot','helt','roll'])
                        mle_recognizers.append(self.mle_shape(i, placements, shape_type))
                    # if it's too short, leave the PSSM unchanged
                    else:
                        mle_recognizers.append(copy.deepcopy(organism.recognizers[i]))
                
                # (III) Leave the PSSM unchanged
                # [with probability 1 - (probability_mle_pssm + probability_mle_pssm_to_shape)]
                else:
                    mle_recognizers.append(copy.deepcopy(organism.recognizers[i]))
            
            # SHAPE
            elif organism.recognizers[i].is_shape():
                x = random.random()
                
                # (I) MLE-optimized Shape
                # [with probability `probability_mle_shape`]
                if x < self.probability_mle_shape:
                    mle_recognizers.append(self.mle_shape(i, placements, organism.recognizers[i].type))
                
                # (II) Turn the Shape into a PSSM (MLE-optimized)
                # [with probability `probability_mle_shape_to_pssm`]
                elif (x - self.probability_mle_shape) < self.probability_mle_shape_to_pssm:
                    mle_recognizers.append(self.mle_pssm(i, placements))
                
                # (III) Leave the Shape unchanged
                # [with probability 1 - (probability_mle_shape + probability_mle_shape_to_pssm)]
                else:
                    mle_recognizers.append(copy.deepcopy(organism.recognizers[i]))
        
        # Return the new organism
        org = self.clone_organism(organism)
        org.recognizers = mle_recognizers
        org.connectors = mle_connectors
        org.flatten()
        return org









