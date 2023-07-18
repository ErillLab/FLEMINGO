from . import constants
import numpy
import itertools
import numpy as np
import math
import pickle






def print_histogram(histogram):
    counter = 0
    print(histogram[1])
    for i,j in zip(histogram[0], histogram[1]):
        counter += 1 
        print('[' + str(round(j, 2)) + ', ' + str(round(histogram[1][counter], 2)) + ') ' + '{0: >10}'.format(i))

def get_null_mgw(sequences, n, num_bins):
    num_pentamers = n - 4
    pentamer_scores = []
    scores = []

    for sequence in sequences:
        
        pentamer_scores = []
        for i in range(num_pentamers):
            index = 0 

            for j in range(i, i + 5):
                
                nuc = sequence[j]
                
                if nuc == 'A':
                    index += pow(4, 4 - (j - i)) * 0
                if nuc == 'G':
                    index += pow(4, 4 - (j - i)) * 1
                if nuc == 'C':
                    index += pow(4, 4 - (j - i)) * 2
                if nuc == 'T':
                    index += pow(4, 4 - (j - i)) * 3

            pentamer_scores.append(constants.MGW_SCORES[index])

        scores.append(sum(pentamer_scores)/num_pentamers)
    
    # Compute frequency for each bin
    counts, edges = np.histogram(scores, bins=num_bins)
    frequencies = counts / sum(counts)
    for i in range(len(frequencies)):
        if frequencies[i] != 0.00:
            frequencies[i] = np.log(frequencies[i])
        else:
            frequencies[i] = np.nan
    
    return (frequencies, edges)

def get_null_prot(sequences, n, num_bins):
    num_pentamers = n - 4
    pentamer_scores = []
    scores = []

    for sequence in sequences:
        
        pentamer_scores = []
        for i in range(num_pentamers):
            index = 0 

            for j in range(i, i + 5):
                
                nuc = sequence[j]
                
                if nuc == 'A':
                    index += pow(4, 4 - (j - i)) * 0
                if nuc == 'G':
                    index += pow(4, 4 - (j - i)) * 1
                if nuc == 'C':
                    index += pow(4, 4 - (j - i)) * 2
                if nuc == 'T':
                    index += pow(4, 4 - (j - i)) * 3

            pentamer_scores.append(constants.PROT_SCORES[index])

        scores.append(sum(pentamer_scores)/num_pentamers)
    
    # Compute frequency for each bin
    counts, edges = np.histogram(scores, bins=num_bins)
    frequencies = counts / sum(counts)
    for i in range(len(frequencies)):
        if frequencies[i] != 0.00:
            frequencies[i] = np.log(frequencies[i])
        else:
            frequencies[i] = np.nan
    
    return (frequencies, edges)

def get_null_roll(sequences, n, num_bins):
    num_pentamers = n - 4
    pentamer_scores = []
    scores = []

    for sequence in sequences:
        
        pentamer_scores = []
        for i in range(num_pentamers):
            index = 0 

            for j in range(i, i + 5):
                
                nuc = sequence[j]
                
                if nuc == 'A':
                    index += pow(4, 4 - (j - i)) * 0
                if nuc == 'G':
                    index += pow(4, 4 - (j - i)) * 1
                if nuc == 'C':
                    index += pow(4, 4 - (j - i)) * 2
                if nuc == 'T':
                    index += pow(4, 4 - (j - i)) * 3

            pentamer_scores.append(constants.ROLL_SCORES[index])
            pentamer_scores.append(constants.ROLL_SCORES[1024 + index])
        
        # Weighted average where first and last element have double weight (counted twice)
        # Note that the number of elements in pentamer_scores is 2*num_pentamers
        scores.append( (pentamer_scores[0] + sum(pentamer_scores) + pentamer_scores[-1]) / ((2*num_pentamers) + 2) )
    
    # Compute frequency for each bin
    counts, edges = np.histogram(scores, bins=num_bins)
    frequencies = counts / sum(counts)
    for i in range(len(frequencies)):
        if frequencies[i] != 0.00:
            frequencies[i] = np.log(frequencies[i])
        else:
            frequencies[i] = np.nan
    
    return (frequencies, edges)

def get_null_helt(sequences, n, num_bins):
    num_pentamers = n - 4
    pentamer_scores = []
    scores = []

    for sequence in sequences:
        
        pentamer_scores = []
        for i in range(num_pentamers):
            index = 0 

            for j in range(i, i + 5):
                
                nuc = sequence[j]
                
                if nuc == 'A':
                    index += pow(4, 4 - (j - i)) * 0
                if nuc == 'G':
                    index += pow(4, 4 - (j - i)) * 1
                if nuc == 'C':
                    index += pow(4, 4 - (j - i)) * 2
                if nuc == 'T':
                    index += pow(4, 4 - (j - i)) * 3

            pentamer_scores.append(constants.HELT_SCORES[index])
            pentamer_scores.append(constants.HELT_SCORES[1024 + index])
        
        # Weighted average where first and last element have double weight (counted twice)
        # Note that the number of elements in pentamer_scores is 2*num_pentamers
        scores.append( (pentamer_scores[0] + sum(pentamer_scores) + pentamer_scores[-1]) / ((2*num_pentamers) + 2) )
    
    # Compute frequency for each bin
    counts, edges = np.histogram(scores, bins=num_bins)
    frequencies = counts / sum(counts)
    for i in range(len(frequencies)):
        if frequencies[i] != 0.00:
            frequencies[i] = np.log(frequencies[i])
        else:
            frequencies[i] = np.nan
    
    return (frequencies, edges)

def generate_range(a, b, models, bins):
    if models == {} or (bins not in models.keys()):
        models[bins] = {}
        models[bins]["mgw"] = {}
        models[bins]["prot"] = {}
        models[bins]["roll"] = {}
        models[bins]["helt"] = {}
    for i in range(a, b):
        print("Generating cartesian product for length {}...".format(i))
        sequences = (list(itertools.product(['A', 'G', 'C', 'T'], repeat=i)))
        print("Finished generating sequences for length {}".format(i))
        print("Calculating null distributions for length {}...".format(i))
        model, edges = get_null_mgw(sequences, i, bins)
        models[bins]["mgw"][i] = {}
        models[bins]["mgw"][i]["frequencies"] = model
        models[bins]["mgw"][i]["bins"] = edges
        model, edges = get_null_prot(sequences, i, bins)
        models[bins]["prot"][i] = {}
        models[bins]["prot"][i]["frequencies"] = model
        models[bins]["prot"][i]["bins"] = edges
        model, edges = get_null_roll(sequences, i, bins)
        models[bins]["roll"][i] = {}
        models[bins]["roll"][i]["frequencies"] = model
        models[bins]["roll"][i]["bins"] = edges
        model, edges = get_null_helt(sequences, i, bins)
        models[bins]["helt"][i] = {}
        models[bins]["helt"][i]["frequencies"] = model
        models[bins]["helt"][i]["bins"] = edges
        print("Finsihed calculating null for shapes of length {}".format(i))
        with open("models/models", "wb") as outfile:
            pickle.dump(models, outfile)
