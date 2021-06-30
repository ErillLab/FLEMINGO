# -*- coding: utf-8 -*-
"""Tests an organism Fitness
"""

from search_organisms import read_fasta_file, read_json_file, export_organism
from objects.organism_factory import OrganismFactory
import json
import numpy as np
import matplotlib.pyplot as plt
import time

CONFIG_FILE = "config.json"

def main():
    """Main execution for the test organisms

    """
    #read configuration file
    config = read_json_file(CONFIG_FILE)
    positive_path = (
        config["main"]["DATASET_BASE_PATH_DIR"]
        + config["main"]["POSITIVE_FILENAME"]
    )
    negative_path = (
        config["main"]["DATASET_BASE_PATH_DIR"]
        + config["main"]["NEGATIVE_FILENAME"]
    )
    max_sequences_to_fit_pos = config["main"]["MAX_SEQUENCES_TO_FIT_POS"]
    max_sequences_to_fit_neg = config["main"]["MAX_SEQUENCES_TO_FIT_NEG"]

    input_organisms_path = config["main"]["INPUT_FILENAME"]
    positive_dataset = read_fasta_file(positive_path)
    #positive_dataset.sort()
    negative_dataset = read_fasta_file(negative_path)
    #print("{} {}".format(len(positive_dataset), len(negative_dataset)))
    
    genome_length = config["main"]["GENOME_LENGTH"]
    
    organism_factory = OrganismFactory(
        config["organism"],
        config["organismFactory"],
        config["connector"],
        config["pssm"],
    )
    
    # create organisms in input file
    a_organisms = organism_factory.import_organisms(input_organisms_path)

    # for each organism: compute Boltzmann, Gini and discrmiminative fitness
    # write out on console the results
    # export the placement of the organism against each pos/neg set sequence
    for org in a_organisms:
        
        # start_time = time.time()

        nodes = org.count_nodes()
        
        # Boltzmannian fitness
        performance1 = org.get_boltz_fitness(positive_dataset[:max_sequences_to_fit_pos],
                                             negative_dataset[:max_sequences_to_fit_neg],
                                             genome_length, use_gini=True)
        boltz_fitness = performance1["score"]

        # Gini coefficient
        gini_coeff = performance1["avg_gini"]
        
        # Kolmogorov fitness
        performance1 = org.get_kolmogorov_fitness(positive_dataset[:max_sequences_to_fit_pos],
                                             negative_dataset[:max_sequences_to_fit_neg],
                                             use_gini=False)
        kolm_fitness = performance1["score"]
        
        
        # print("Boltzmann --- %s seconds ---" % (time.time() - start_time))
        
        # start_time = time.time()
        
        # Discriminative fitness
        Pf = org.get_additive_fitness(positive_dataset[:max_sequences_to_fit_pos],
                                      use_gini=False)
        P = Pf["score"]
        Pstd = Pf["stdev"]
        
        Nf = org.get_additive_fitness(negative_dataset[:max_sequences_to_fit_neg],
                                      use_gini=False)
        N = Nf["score"]
        Nstd = Nf["stdev"]
        
        discr_fitness =  P - N
        # print("Additive --- %s seconds ---" % (time.time() - start_time))
        
        
        # Welch's fitness
        if Pstd < 1:
            Pstd_for_welchs = 1
        else:
            Pstd_for_welchs = Pstd
        
        if Nstd < 1:
            Nstd_for_welchs = 1
        else:
            Nstd_for_welchs = Nstd
        
        sterr_P = Pstd_for_welchs / max_sequences_to_fit_pos**(1/2)
        sterr_N = Nstd_for_welchs / max_sequences_to_fit_neg**(1/2)
        welchs_fitness =  (P - N) / (sterr_P**2 + sterr_N**2)**(1/2)
        
        
        # Distribution plotting
        Pe = org.get_binding_energies(positive_dataset[:max_sequences_to_fit_pos],
                                     traceback=False, use_gini=False)
        
        Ne = org.get_binding_energies(negative_dataset[:max_sequences_to_fit_neg],
                                      traceback=False, use_gini=False)
        
        plt.hist(Pe, alpha=0.5, label='positive set')
        plt.hist(Ne, alpha=0.5, label='negative set')
        plt.legend()
        plt.savefig('org_' + str(org._id) + '_energy_distr')
        plt.close()
        
        # print out results (fitness, nodes, Gini, etc.) for organism
        print(
            (
                "Org {} Nodes: {:.2f} GiniPSSMs: {:.2f} P: {:.2f}(+/-){:.2f} N: {:.2f}(+/-){:.2f}"
                + " DiscrF: {:.2f} BoltzF: {:.2f} KolmF: {:.2f} WelchsF: {:.2f}\n"
            ).format(
                org._id,  # "Org"
                nodes,  # "Nodes"
                gini_coeff, # GiniPSSMs
                P,  # "P"
                Pstd, # stdev for P
                N,  # "N"
                Nstd, # stdev for N
                discr_fitness,  # "DiscrF"
                boltz_fitness,  # "BoltzF"
                kolm_fitness,    # "KolmF"
                welchs_fitness  # "WelchF"
                )
        )
        
                
        #export the organism results on the positive dataset
        export_organism(
            org,
            positive_dataset,
            "{}positive_{}".format(
                config["main"]["RESULT_TEST_BASE_PATH_DIR"], org._id
            ),
            organism_factory,
        )

        #export the organism results on the negative dataset
        export_organism(
            org,
            negative_dataset[:50],
            "{}negative_{}".format(
                config["main"]["RESULT_TEST_BASE_PATH_DIR"], org._id
            ),
            organism_factory,
        )
    
    print("\n(Energy histograms are saved in the current working directory)")
        

if __name__ == "__main__":
    
    # TODO: Add the profiler here to improve the fitness calculation performance
    
    main()



