# -*- coding: utf-8 -*-

"""Main execution

This program searches for models that fit a specific motif.

"""

import time
import random
import copy
import json
import os
import cProfile
import pstats
import io
import numpy as np
from objects.organism_factory import OrganismFactory
from Bio import SeqIO

"""
Variable definition
"""
POPULATION_LENGTH = 0
DATASET_BASE_PATH_DIR = ""
RESULT_BASE_PATH_DIR = ""
POSITIVE_FILENAME = ""
NEGATIVE_FILENAME = ""
POPULATION_ORIGIN = ""
POPULATION_FILL_TYPE = ""
INPUT_FILENAME = ""
OUTPUT_FILENAME = ""
MAX_SEQUENCES_TO_FIT_POS = 0
MAX_SEQUENCES_TO_FIT_NEG = 0
MIN_ITERATIONS = 0
MIN_FITNESS = 0
RECOMBINATION_PROBABILITY = 0.0
THRESHOLD = 0.0

JSON_CONFIG_FILENAME = "config.json"
"""
Configuration of the object types
Populated by JSON read.
"""
configOrganism: dict = {}
configOrganismFactory: dict = {}
configConnector: dict = {}
configPssm: dict = {}


# list that holds the population
organism_population: list = []

# mean_nodes are the average number of nodes per organism in the population
# used to calculate organism complexity
mean_nodes: float = 0
# mean_fitness is the average fitness per organism in the population
# used to calculate organism complexity
mean_fitness: float = 0

# Initialize datasets
positive_dataset: list = []
negative_dataset: list = []


def main():
    """Main function for the motif seek
    """

    print("Loading parameters...")
    positive_dataset = read_fasta_file(
        DATASET_BASE_PATH_DIR + POSITIVE_FILENAME
    )
    negative_dataset = read_fasta_file(
        DATASET_BASE_PATH_DIR + NEGATIVE_FILENAME
    )

    mean_nodes = 0
    mean_fitness = 0

    print("Instantiating population...")

    """
    Generate initial population
    """
    # Instantiate organism Factory object with object configurations
    organism_factory = OrganismFactory(
        configOrganism, configOrganismFactory, configConnector, configPssm
    )
    # initialize organisms population
    organism_population = []

    # Generate population depending on origin and fill type.
    # Origin can be "random" or "file"(read a set of organisms from a file).
    if POPULATION_ORIGIN.lower() == "random":
        # For a random origin, we can generate #POPULATION_LENGTH organisms.
        for i in range(POPULATION_LENGTH):
            new_organism = organism_factory.get_organism()
            organism_population.append(new_organism)
    elif POPULATION_ORIGIN.lower() == "file":
        #if the population is seeded from file
        # Set the file organisms and fill with random/same organisms
        # POPULATION_LENGTH must be >= len(fileOrganisms)
        file_organisms = organism_factory.import_organisms(INPUT_FILENAME)
        remaining_organisms = POPULATION_LENGTH - len(file_organisms)
        fill_organism_population = []

        if POPULATION_FILL_TYPE.lower() == "random":
            # FILL remainder of the population WITH RANDOM organisms
            for i in range(remaining_organisms):
                new_organism = organism_factory.get_organism()
                fill_organism_population.append(new_organism)

        elif POPULATION_FILL_TYPE.lower() == "uniform":
            # FILL the remainder of the population WITH ORGANISMS IN FILE
            for i in range(remaining_organisms):
                new_organism = copy.deepcopy(
                    file_organisms[i % len(file_organisms)]
                )
                fill_organism_population.append(new_organism)
                new_organism.set_id(organism_factory.get_id())

        # join 
        organism_population = file_organisms + fill_organism_population

    else:
        raise Exception("Not a valid population origin, "
            + "check the configuration file.")
    

    print("Population size = {}".format(len(organism_population)))
    """
    Initialize iteration variables.
    """
    iterations = 0
    max_score = float("-inf")
    last_max_score = 0.0
    # Organism with highest fitness in the simulation
    best_organism = (
        None,  # the organism object
        -np.inf,  # its fitness
        0,  # number of nodes it is made of
        0.0  # the penalty it gets
    )
    # Organism with highest fitness in the iteration
    max_organism = (
        None,  # the organism object
        -np.inf,  # its fitness
        0,  # number of nodes it is made of
        0.0  # the penalty it gets
    )
    timeformat = "%Y-%m-%d--%H-%M-%S"
    print("Starting execution...")

    # Main loop, it iterates until organisms do not get a significant change
    # or MIN_ITERATIONS or MIN_FITNESS is reached.

    while not is_finished(END_WHILE_METHOD, iterations, max_score, 
                          last_max_score):
        
        # Shuffle population
        # Organisms are shuffled for deterministic crowding selection
        random.shuffle(organism_population)
        
        # Shuffle datasets
        # Datasets are shuffled for subsampling
        if RANDOM_SHUFFLE_SAMPLING_POS:
            random.shuffle(positive_dataset)
        if RANDOM_SHUFFLE_SAMPLING_NEG:
            random.shuffle(negative_dataset)

        # Reset max_score
        last_max_score = max_score
        max_score = float("-inf")
        changed_best_score = False
        initial = time.time()
        
        a_fitness = []
        a_nodes = []

        # Deterministic crowding
        # Iterate over pairs of organisms
        for i in range(0, len(organism_population) - 1, 2):
            org1 = organism_population[i]
            org2 = organism_population[i + 1]
            
            pos_set_sample = random.sample(positive_dataset, 3)  # !!! Temporarily hardcoded number of sequences
            ref_seq = random.choice(pos_set_sample)
            
            # Cross parents to get children
            # Recombination process
            child1, child2 = organism_factory.get_children(
                org1, org2, ref_seq, pos_set_sample
            )
            
            # Make two pairs: each parent is paired with the more similar child
            # (the child with higher ratio of nodes from that parent).
            pair_children = []
            ''' pair_children is a list of two elements. The two parents we are
            now working with are elements i and i+1 in  organism_population.
            We need to pair each of them with one of the two children obtained
            with  get_children  method.
            
                - The first element in  pair_children  will be a tuple where
                  the first element is organism i, and the second one is a
                  child
                  
                - The second element in  pair_children  will be a tuple where
                  the first element is organism i+1, and the second one is the
                  other child
            '''
            
            # get parent1/parent2 ratio for the children
            child1_p1p2_ratio = child1.get_parent1_parent2_ratio()
            child2_p1p2_ratio = child2.get_parent1_parent2_ratio()
            
            # If a parent gets paired with an empty child, the empty child is
            # substituted by a deepcopy of the parent, i.e. the parent escapes
            # competition
            if child1_p1p2_ratio > child2_p1p2_ratio:
                # org1 with child1
                if child1.count_nodes() > 0:
                    pair_children.append( (org1, child1) )
                else:
                    pair_children.append( (org1, copy.deepcopy(org1)) )
                # org2 with child2
                if child2.count_nodes() > 0:
                    pair_children.append( (org2, child2) )
                else:
                    pair_children.append( (org2, copy.deepcopy(org2)) )
            else:
                # org1 with child2
                if child2.count_nodes() > 0:
                    pair_children.append( (org1, child2) )
                else:
                    pair_children.append( (org1, copy.deepcopy(org1)) )
                # org2 with child1
                if child1.count_nodes() > 0:
                    pair_children.append( (org2, child1) )
                else:
                    pair_children.append( (org2, copy.deepcopy(org2)) )
            
            # Mutate the children
            child1.mutate(organism_factory)
            child2.mutate(organism_factory)
            
            # Make the two organisms in each pair compete
            # j index is used to re insert winning organism into the population
            for j in range(len(pair_children)):
                '''
                when j is 0 we are dealing with organism i
                when j is 1 we are dealing with organism i+1
                
                Therefore, the winner of the competition will replace element
                i+j in  organism_population
                '''

                first_organism = pair_children[j][0]  # Parent Organism
                second_organism = pair_children[j][1]  # Child Organism
                
                # Boltzmannian fitness
                if FITNESS_FUNCTION == "boltzmannian":
                    performance1 = first_organism.get_boltz_fitness(positive_dataset[:MAX_SEQUENCES_TO_FIT_POS],
                                                                    negative_dataset[:MAX_SEQUENCES_TO_FIT_NEG],
                                                                    GENOME_LENGTH, use_gini=False)
                    fitness1 = performance1["score"]
                    gini1 = performance1["avg_gini"]
                    
                    performance2 = second_organism.get_boltz_fitness(positive_dataset[:MAX_SEQUENCES_TO_FIT_POS],
                                                                     negative_dataset[:MAX_SEQUENCES_TO_FIT_NEG],
                                                                     GENOME_LENGTH, use_gini=False)
                    fitness2 = performance2["score"]
                    gini2 = performance2["avg_gini"]
                    
                    fitness1 = round(fitness1, 8)
                    fitness2 = round(fitness2, 8)
                
                # Kolmogorov fitness
                # Computes Kolmogorov-Smirnov test on positive/negative set scores
                elif FITNESS_FUNCTION == "kolmogorov":
                    performance1 = first_organism.get_kolmogorov_fitness(positive_dataset[:MAX_SEQUENCES_TO_FIT_POS],
                                                                    negative_dataset[:MAX_SEQUENCES_TO_FIT_NEG],
                                                                    use_gini=False)
                    fitness1 = performance1["score"]
                    gini1 = performance1["avg_gini"]                   
                    
                    performance2 = second_organism.get_kolmogorov_fitness(positive_dataset[:MAX_SEQUENCES_TO_FIT_POS],
                                                                    negative_dataset[:MAX_SEQUENCES_TO_FIT_NEG],
                                                                    use_gini=False)

                    fitness2 = performance2["score"]
                    gini2 = performance2["avg_gini"]
                    
                    fitness1 = round(fitness1, 8)
                    fitness2 = round(fitness2, 8)

                # Discriminative fitness
                elif FITNESS_FUNCTION == "discriminative":
                    positive_performance1 = first_organism.get_additive_fitness(positive_dataset[:MAX_SEQUENCES_TO_FIT_POS],
                                                                        use_gini=False)
                    negative_performance1 = first_organism.get_additive_fitness(negative_dataset[:MAX_SEQUENCES_TO_FIT_NEG],
                                                                        use_gini=False)
                    p_1 = positive_performance1["score"]
                    n_1 = negative_performance1["score"]
                    fitness1 =  p_1 - n_1
                    gini1 = positive_performance1["avg_gini"]
                    
                    positive_performance2 = second_organism.get_additive_fitness(positive_dataset[:MAX_SEQUENCES_TO_FIT_POS],
                                                                        use_gini=False)
                    negative_performance2 = second_organism.get_additive_fitness(negative_dataset[:MAX_SEQUENCES_TO_FIT_NEG],
                                                                        use_gini=False)
                    p_2 = positive_performance2["score"]
                    n_2 = negative_performance2["score"]
                    fitness2 =  p_2 - n_2
                    gini2 = positive_performance2["avg_gini"]
                
                elif FITNESS_FUNCTION == "welchs":
                    # First organism
                    positive_performance1 = first_organism.get_additive_fitness(positive_dataset[:MAX_SEQUENCES_TO_FIT_POS],
                                                                        use_gini=False)
                    negative_performance1 = first_organism.get_additive_fitness(negative_dataset[:MAX_SEQUENCES_TO_FIT_NEG],
                                                                        use_gini=False)
                    p_1 = positive_performance1["score"]
                    n_1 = negative_performance1["score"]
                    
                    # Standard deviations
                    sigma_p_1 = positive_performance1["stdev"]
                    sigma_n_1 = negative_performance1["stdev"]
                    
                    # Lower bound to sigma
                    # (Being more consistent than that on the sets will not help
                    # your fitness)
                    if sigma_p_1 < 1:
                        sigma_p_1 = 1
                    if sigma_n_1 < 1:
                        sigma_n_1 = 1
                    
                    # Standard errors
                    sterr_p_1 = sigma_p_1 / MAX_SEQUENCES_TO_FIT_POS**(1/2)
                    sterr_n_1 = sigma_n_1 / MAX_SEQUENCES_TO_FIT_NEG**(1/2)
                    
                    # Welch's t score
                    fitness1 =  (p_1 - n_1) / (sterr_p_1**2 + sterr_n_1**2)**(1/2)
                    
                    gini1 = positive_performance1["avg_gini"]
                    
                    # Second organism
                    positive_performance2 = second_organism.get_additive_fitness(positive_dataset[:MAX_SEQUENCES_TO_FIT_POS],
                                                                        use_gini=False)
                    negative_performance2 = second_organism.get_additive_fitness(negative_dataset[:MAX_SEQUENCES_TO_FIT_NEG],
                                                                        use_gini=False)
                    p_2 = positive_performance2["score"]
                    n_2 = negative_performance2["score"]
                    
                    # Standard deviations
                    sigma_p_2 = positive_performance2["stdev"]
                    sigma_n_2 = negative_performance2["stdev"]
                    
                    # Lower bound to sigma
                    # (Being more consistent than that on the sets will not help
                    # your fitness)
                    if sigma_p_2 < 1:
                        sigma_p_2 = 1
                    if sigma_n_2 < 1:
                        sigma_n_2 = 1
                    
                    # Standard errors
                    sterr_p_2 = sigma_p_2 / MAX_SEQUENCES_TO_FIT_POS**(1/2)
                    sterr_n_2 = sigma_n_2 / MAX_SEQUENCES_TO_FIT_NEG**(1/2)
                    
                    # Welch's t score
                    fitness2 =  (p_2 - n_2) / (sterr_p_2**2 + sterr_n_2**2)**(1/2)
                    
                    gini2 = positive_performance2["avg_gini"]
                
                else:
                    raise Exception("Not a valid fitness function name, "
                                    + "check the configuration file.")


                if MAX_NODES != None:  # Upper_bound to complexity
                    
                    if first_organism.count_nodes() > MAX_NODES:
                        #print(first_organism.count_nodes(), "nodes")
                        fitness1 = -1000 * int(first_organism.count_nodes())
                    
                    if second_organism.count_nodes() > MAX_NODES:
                        #print(second_organism.count_nodes(), "nodes")
                        fitness2 = -1000 * int(second_organism.count_nodes())

                if MIN_NODES != None:  # Lower_bound to complexity
                    
                    if first_organism.count_nodes() < MIN_NODES:
                        #print(first_organism.count_nodes(), "nodes")
                        fitness1 = -1000 * int(first_organism.count_nodes())
                    
                    if second_organism.count_nodes() < MIN_NODES:
                        #print(second_organism.count_nodes(), "nodes")
                        fitness2 = -1000 * int(second_organism.count_nodes())                
                
                if INEQUALITY_PENALTY_METHOD=="avg_gini":
                    # INEQUALITY_PENALTY_PARAM acts as a penalty buffer
                    # The higher this parameter, the less important is the effect of the Gini penalty
                    # It's meant to work in the range [1, +inf)
                    effective_fitness_1 = fitness1 * (INEQUALITY_PENALTY_PARAM - gini1)
                    effective_fitness_2 = fitness2 * (INEQUALITY_PENALTY_PARAM - gini2)
                else:
                    effective_fitness_1 = fitness1
                    effective_fitness_2 = fitness2
                
                
                if (
                        effective_fitness_1 > effective_fitness_2
                ):  # The first organism wins
                    # Set it back to the population and save fitness
                    # for next iteration
                    organism_population[i + j] = first_organism
                    a_fitness.append(fitness1)
                    # If the parent wins, mean_nodes don't change
                    a_nodes.append(first_organism.count_nodes())

                    # Check if its the max score in that iteration
                    if effective_fitness_1 > max_score:
                        max_score = effective_fitness_1
                        max_organism = (
                            first_organism,
                            effective_fitness_1,
                            first_organism.count_nodes(),
                            gini1,
                        )

                    # Check if its the max score so far and if it is set it as
                    # best organism
                    if max_organism[1] > best_organism[1]:
                        # ID, EF, Nodes, Penalty applied
                        best_organism = max_organism
                        changed_best_score = True

                else:  # The second organism wins (child)
                    # Set it back to the population and save fitness for next
                    # iteration
                    organism_population[i + j] = second_organism
                    a_fitness.append(fitness2)
                    # If the child wins, update mean_nodes
                    # mean_nodes = ((meanNodes * POPULATION_LENGTH) +
                    # second_organism.count_nodes() -
                    # first_organism.count_nodes()) / POPULATION_LENGTH
                    a_nodes.append(second_organism.count_nodes())

                    # Check if its the max score in that iteration
                    if effective_fitness_2 > max_score:
                        max_score = effective_fitness_2
                        max_organism = (
                            second_organism,
                            effective_fitness_2,
                            second_organism.count_nodes(),
                            gini2,
                        )

                    # Check if its the max score so far and if it is set it as
                    # best organism
                    if effective_fitness_2 > best_organism[1]:
                        # ID, EF, Nodes, Penalty applied
                        best_organism = max_organism
                        changed_best_score = True

                # END FOR j

            # END FOR i

        # Mean fitness in the population
        mean_fitness = np.mean(a_fitness)
        # Standard deviation of fitness in the population
        standard_dev_fitness = np.std(a_fitness)
        # Inequality of fitness in the population (measured with the Gini coefficient)
        gini_fitness = gini_RSV(a_fitness)
        # Mean number of nodes per organism in the population
        mean_nodes = np.mean(a_nodes)

        # Show IDs of final array
        # print("-"*10)
        _m, _s = divmod((time.time() - initial), 60)
        _h, _m = divmod(_m, 60)
        s_time = "{}h:{}m:{:.2f}s".format(int(_h), int(_m), _s)
        print_ln(
            (
                "Iter: {} AF:{:.2f} SDF:{:.2f} GF:{:.2f} AN:{:.2f}"
                + " - MO: {} MF: {:.2f} MN: {} MP: {:.2f}"
                + " -  BO: {} BF: {:.2f} BN: {} BP: {:.2f} Time: {}"
            ).format(
                iterations,  # "Iter"
                mean_fitness,  # "AF"
                standard_dev_fitness,  # "SDF"
                gini_fitness,  # "GF"
                mean_nodes,  # "AN"
                max_organism[0]._id,  # "MO"
                max_organism[1],  # "MF" (fitness)
                max_organism[2],  # "MN" (nodes)
                max_organism[3],  # "MP" (penalty)
                best_organism[0]._id,  # "BO"
                best_organism[1],  # "BF" (fitness)
                best_organism[2],  # "BN" (nodes)
                best_organism[3],  # "BP" (penalty)
                s_time,  # Time
            ),
            RESULT_BASE_PATH_DIR + OUTPUT_FILENAME,
        )
        
        # Print against a random positive sequence
        pos_seq_index = random.randint(0, len(positive_dataset)-1)
        placement = max_organism[0].get_placement(positive_dataset[pos_seq_index], traceback=True)
        placement.print_placement(stdout = True)
        
        if RANDOM_SHUFFLE_SAMPLING_POS:
            # Sort the positive dataset so that, when exporting with  export_organism
            # function, the subset of DNA sequences used for printing will be
            # always the same (and they'll apear always in the same order).
            positive_dataset.sort()
        
        # Export organism if new best organism
        if changed_best_score:
            filename = "{}_{}".format(
                time.strftime(timeformat), best_organism[0]._id
            )
            export_organism(
                best_organism[0], positive_dataset, filename, organism_factory
            )
        # Periodic organism export
        if iterations % PERIODIC_EXPORT == 0:
            filename = "{}_{}".format(
                time.strftime(timeformat), max_organism[0]._id
            )
            export_organism(
                max_organism[0], positive_dataset, filename, organism_factory
            )
        
        iterations += 1
        # END WHILE

    # TODO: Maybe a good idea to export the full population after all
    # organism_factory.export_organisms(organism_population,
    #         RESULT_BASE_PATH_DIR+"final_population.json")


def is_finished(
        method: str, iterations: int, max_score: float, last_max_score: float
) -> bool:
    """Checks if main while loop is finished
    methods: 'Iterations', 'minScore', 'Threshold'

    Args:
        method: Name of the finishing method
        max_score: max score recorded on the current iteration
        last_max_score: max score recorded on the laset iteration
        iterations: Number of the current iteration

    Returns:
        True if program should finnish.
        False otherwise
    """

    if method.lower() == "iterations":
        return iterations >= MIN_ITERATIONS

    if method.lower() == "fitness":
        return max_score >= MIN_FITNESS

    if method.lower() == "threshold":
        return abs(last_max_score - max_score) <= THRESHOLD

    return True


def export_organism(
        organism, dataset: list, filename: str, factory: OrganismFactory
) -> None:
    """Exports a single organism in json format, visual format and its
    recognizers binding

    Args:
        organism (OrganismObject): organism to export
        dataset: Sequences to check the organism binding
        filename: Previous info to export filenames. Common in all filenames
        factory: Used to export in json format
    """

    organism_file = "{}{}_organism.txt".format(RESULT_BASE_PATH_DIR, filename)
    organism_file_json = "{}{}_organism.json".format(
        RESULT_BASE_PATH_DIR, filename
    )
    results_file = "{}{}_results.txt".format(RESULT_BASE_PATH_DIR, filename)

    organism.export(organism_file)
    organism.export_results(dataset, results_file)
    factory.export_organisms([organism], organism_file_json)


def set_up():
    """Reads configuration file and sets up all program variables

    """

    # specify as global variable so it can be accesed in local
    # contexts outside setUp

    global END_WHILE_METHOD
    global POPULATION_LENGTH
    global DATASET_BASE_PATH_DIR
    global RESULT_BASE_PATH_DIR
    global POSITIVE_FILENAME
    global NEGATIVE_FILENAME
    global RESULT_PATH_PATH_DIR
    global MAX_SEQUENCES_TO_FIT_POS
    global MAX_SEQUENCES_TO_FIT_NEG
    global RANDOM_SHUFFLE_SAMPLING_POS
    global RANDOM_SHUFFLE_SAMPLING_NEG
    global FITNESS_FUNCTION
    global USE_GINI
    global GENOME_LENGTH
    global INEQUALITY_PENALTY_METHOD
    global INEQUALITY_PENALTY_PARAM
    global MIN_ITERATIONS
    global MIN_FITNESS
    global THRESHOLD
    global POPULATION_ORIGIN
    global POPULATION_FILL_TYPE
    global INPUT_FILENAME
    global OUTPUT_FILENAME
    global PERIODIC_EXPORT
    global MAX_NODES
    global MIN_NODES

    # Config data
    global configOrganism
    global configOrganismFactory
    global configConnector
    global configPssm

    config = read_json_file(JSON_CONFIG_FILENAME)
    
    # Store config variables for main function
    
    POPULATION_LENGTH = config["main"]["POPULATION_LENGTH"]
    if POPULATION_LENGTH % 2 == 1:
        raise Exception("POPULATION_LENGTH must be an even number. It's " +
                        "currently set to " + str(POPULATION_LENGTH))
    DATASET_BASE_PATH_DIR = config["main"]["DATASET_BASE_PATH_DIR"]
    RESULT_BASE_PATH_DIR = (
        config["main"]["RESULT_BASE_PATH_DIR"]
        + time.strftime("%Y%m%d%H%M%S")
        + "/"
    )
    POSITIVE_FILENAME = config["main"]["POSITIVE_FILENAME"]
    NEGATIVE_FILENAME = config["main"]["NEGATIVE_FILENAME"]
    MAX_SEQUENCES_TO_FIT_POS = config["main"]["MAX_SEQUENCES_TO_FIT_POS"]
    MAX_SEQUENCES_TO_FIT_NEG = config["main"]["MAX_SEQUENCES_TO_FIT_NEG"]
    RANDOM_SHUFFLE_SAMPLING_POS = config["main"]["RANDOM_SHUFFLE_SAMPLING_POS"]
    RANDOM_SHUFFLE_SAMPLING_NEG = config["main"]["RANDOM_SHUFFLE_SAMPLING_NEG"]
    FITNESS_FUNCTION = config["main"]["FITNESS_FUNCTION"]
    USE_GINI = config["main"]["USE_GINI"]
    GENOME_LENGTH = config["main"]["GENOME_LENGTH"]
    INEQUALITY_PENALTY_METHOD = config["main"]["INEQUALITY_PENALTY_METHOD"]
    INEQUALITY_PENALTY_PARAM = config["main"]["INEQUALITY_PENALTY_PARAM"]
    MIN_ITERATIONS = config["main"]["MIN_ITERATIONS"]
    MIN_FITNESS = config["main"]["MIN_FITNESS"]
    THRESHOLD = config["main"]["THRESHOLD"]
    END_WHILE_METHOD = config["main"]["END_WHILE_METHOD"]
    POPULATION_ORIGIN = config["main"]["POPULATION_ORIGIN"]
    POPULATION_FILL_TYPE = config["main"]["POPULATION_FILL_TYPE"]
    INPUT_FILENAME = config["main"]["INPUT_FILENAME"]
    OUTPUT_FILENAME = config["main"]["OUTPUT_FILENAME"]
    PERIODIC_EXPORT = config["main"]["PERIODIC_EXPORT"]
    MAX_NODES = config["organism"]["MAX_NODES"]
    MIN_NODES = config["organism"]["MIN_NODES"]

    # Create directory where the output and results will be stored
    os.mkdir(RESULT_BASE_PATH_DIR)

    # Store Config into variables to use later
    configOrganism = config["organism"]
    configOrganismFactory = config["organismFactory"]
    configConnector = config["connector"]
    configPssm = config["pssm"]

    # Throw config on a file
    parameters_path = RESULT_BASE_PATH_DIR + "parameters.txt"
    print_ln("-" * 50, parameters_path)
    print_ln(" " * 20 + "PARAMETERS", parameters_path)
    print_ln("-" * 50, parameters_path)

    print_config_json(config["main"], "Main Config", parameters_path)
    print_config_json(configOrganism, "Organism Config", parameters_path)
    print_config_json(
        configOrganismFactory, "Organism Factory Config", parameters_path
    )
    print_config_json(configConnector, "Connector Config", parameters_path)
    print_config_json(configPssm, "PSSM Config", parameters_path)

    print_ln("-" * 50, parameters_path)


def read_fasta_file(filename: str) -> list:
    """Reads a fasta file and returns an array of DNA sequences (strings)

    TODO: probably it can be useful to create our own Sequence object that
    creates the string and stores some properties from fasta format. Also
    we can adapt the current program to use Biopythons's Seq object.

    Args:
        filename: Name of the file that contains FASTA format sequences to read

    Returns:
        The set of sequences in string format

    """
    dataset = []

    fasta_sequences = SeqIO.parse(open(filename), "fasta")

    for fasta in fasta_sequences:
        dataset.append(str(fasta.seq).lower())

    return dataset


def read_json_file(filename: str) -> dict:
    """Reads a JSON file and returns a dictionary with the content

    Args:
        filename: Name of the json file to read

    Returns:
        Dictionary with the json file info

    """

    with open(filename) as json_content:

        return json.load(json_content)


def print_config_json(config: dict, name: str, path: str) -> None:
    """Print the config file on std out and send it to a file.
    It is useful so we can know which was the configuration on every run

    Args:
        config: Configuration file to print
        name: Title for the configuration file
        path: File to export the configuration info
    """
    print_ln("{}:".format(name), path)

    for key in config.keys():
        print_ln("{}: {}".format(key, config[key]), path)
    print_ln("\n", path)


def print_ln(string: str, name_file: str) -> None:
    """Shows the string on stdout and write it to a file
    (like the python's logging modules does)

    Args:
        string: Information to print on stdout and file
        name_file: path to the file to export the string
    """

    print(string)

    # Here we are sure file exists
    _f = open(name_file, "a+")
    _f.write(string + "\n")
    _f.close()


def gini_RSV(values_for_each_class):
    '''
    Gini coefficient, modified in order to be alble to deal with negative
    values as in "Inequality measures and the issue of negative incomes"
    (Raffinetti, Siletti, Vernizzi)

    Parameters
    ----------
    values_for_each_class : array-like object
        Values associated to each class.
        They don't need to be already sorted and/or normalized.
        They can also be negative.

    Returns
    -------
    giniRSV : float
        Ranges from 0 (perfect equality) to 1 (maximal inequality).

    '''
    
    N = len(values_for_each_class)
    
    numerator = 0
    for i in values_for_each_class:
        for j in values_for_each_class:
            numerator += abs(i - j)
    
    pos = 0  # sum over the positive values
    neg = 0  # sum over the negative values (in absolute value)
    for x in values_for_each_class:
        if x >= 0:
            pos += x
        else:
            neg += -x
    
    mu_RSV = (N - 1) * (pos + neg) / N**2  # modified mu parameter
    
    if mu_RSV == 0:
        # Manage two special cases (avoiding 0-division error):
        #   - when a single value is the input
        #   - when all the values in the input are 0
        # In both cases mu_RSV will be 0
        # No inequality is measurable, and therefore 0 is returned
        return 0
    denominator = 2 * N**2 * mu_RSV
    giniRSV = numerator / denominator
    
    return giniRSV


# Entry point to app execution
# It calculates the time, but could include other app stats

if __name__ == "__main__":

    INITIAL = time.time()
    # Reads configuration file and sets up all program variables
    set_up()

    # Profiling
    PROFILER = cProfile.Profile()
    PROFILER.enable()

    # Main function
    main()

    # Profiler output
    PROFILER.disable()
    STRING = io.StringIO()
    SORTBY = "cumulative"
    PS = pstats.Stats(PROFILER, stream=STRING).sort_stats(SORTBY)
    PS.print_stats()
    # print(STRING.getvalue())

    # Print final execution time and read parameters
    _M, _S = divmod((time.time() - INITIAL), 60)
    _H, _M = divmod(_M, 60)
    print_ln(
        "Time: {}h:{}m:{:.2f}s".format(int(_H), int(_M), _S),
        RESULT_BASE_PATH_DIR + "parameters.txt",
    )
