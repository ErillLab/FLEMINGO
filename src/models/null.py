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
    counter = n

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

        scores.append(sum(pentamer_scores)/len(pentamer_scores))

    hist = np.histogram(scores, bins=num_bins, density=True)
    for i in range(len(hist[0])):
        if hist[0][i] != 0.00:
            hist[0][i] = np.log(hist[0][i])
        else:
            hist[0][i] = np.nan
            
    return hist

def get_null_prot(sequences, n, num_bins):
    num_pentamers = n - 4
    pentamer_scores = []
    scores = []
    counter = n

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

        scores.append(sum(pentamer_scores)/len(pentamer_scores))
    hist = np.histogram(scores, bins=num_bins, density=True)

    for i in range(len(hist[0])):
        if hist[0][i] != 0.00:
            hist[0][i] = np.log(hist[0][i])
        else:
            hist[0][i] = np.nan
    return hist

def get_null_roll(sequences, n, num_bins):
    num_pentamers = n - 4
    pentamer_scores = []
    scores = []
    counter = n

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

        scores.append( (pentamer_scores[0] + sum(pentamer_scores) + pentamer_scores[-1]) / (len(pentamer_scores) + 2) )
    hist = np.histogram(scores, bins=num_bins, density=True)

    for i in range(len(hist[0])):
        if hist[0][i] != 0.00:
            hist[0][i] = np.log(hist[0][i])
        else:
            hist[0][i] = np.nan

    return hist

def get_null_helt(sequences, n, num_bins):
    num_pentamers = n - 4
    pentamer_scores = []
    scores = []
    counter = n

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

        scores.append((pentamer_scores[0] + sum(pentamer_scores) + pentamer_scores[-1]) / (len(pentamer_scores) + 2))
    hist = np.histogram(scores, bins=num_bins, density=True)
    for i in range(len(hist[0])):
        if hist[0][i] != 0.00:
            hist[0][i] = np.log(hist[0][i])
        else:
            hist[0][i] = np.nan
    return hist

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
