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
        print('[' + str(round(j, 1)) + ', ' + str(round(histogram[1][counter], 1)) + ') ' + '{0: >10}'.format(i))

def get_null_mgw(sequences, n):
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
    min_range = round(round(min(scores), 1) - 0.1, 1)
    max_range = round(round(max(scores), 1) + 0.1, 1)
    num_bins = math.ceil(round((max_range - min_range), 1) * 10)
    if (min(scores) > min_range + 0.1):
        num_bins -= 1
        min_range += 0.1

    if (max(scores) < max_range - 0.1):
        num_bins -= 1
        max_range -= 0.1

    return np.histogram(scores, bins=num_bins, density=True, range=(round(min_range, 2), round(max_range, 2)))

def get_null_prot(sequences, n):
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
    min_range = round(round(min(scores), 1) - 0.1, 1)
    max_range = round(round(max(scores), 1) + 0.1, 1)
    num_bins = math.ceil((max_range - min_range) * 10)
    if (min(scores) > min_range + 0.1):
        num_bins -= 1
        min_range += 0.1

    if (max(scores) < max_range - 0.1):
        num_bins -= 1
        max_range -= 0.1
    return np.histogram(scores, bins=num_bins, density=True, range=(min_range, max_range))

def get_null_roll(sequences, n):
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
    min_range = round(round(min(scores), 1) - 0.1, 1)
    max_range = round(round(max(scores), 1) + 0.1, 1)
    num_bins = math.ceil(round((max_range - min_range), 1) * 10)
    if (min(scores) > min_range + 0.1):
        num_bins -= 1
        min_range += 0.1

    if (max(scores) < max_range - 0.1):
        num_bins -= 1
        max_range -= 0.1
    return np.histogram(scores, bins=num_bins, density=True, range=(min_range, max_range))

def get_null_helt(sequences, n):
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

        scores.append( (pentamer_scores[0] + sum(pentamer_scores) + pentamer_scores[-1]) / (len(pentamer_scores) + 2) )
    min_range = round(round(min(scores), 1) - 0.1, 1)
    max_range = round(round(max(scores), 1) + 0.1, 1)
    num_bins = math.ceil((max_range - min_range) * 10)

    if (min(scores) > min_range + 0.1):
        num_bins -= 1
        min_range += 0.1

    if (max(scores) < max_range - 0.1):
        num_bins -= 1
        max_range -= 0.1
    return np.histogram(scores, bins=num_bins, density=True, range=(min_range, max_range))


def generate_range(a, b, models):
    if models == {}:
        models["mgw"] = {}
        models["prot"] = {}
        models["roll"] = {}
        models["helt"] = {}
    for i in range(a, b):
        print("Generating cartesian product for length {}...".format(i))
        sequences = (list(itertools.product(['A', 'G', 'C', 'T'], repeat=i)))
        print("Finished generating sequences for length {}".format(i))
        print("Calculating null distributions for length {}...".format(i))
        model, edges = get_null_mgw(sequences, i)
        models["mgw"][i] = {}
        models["mgw"][i]["frequencies"] = model
        models["mgw"][i]["bins"] = edges
        model, edges = get_null_helt(sequences, i)
        models["prot"][i] = {}
        models["prot"][i]["frequencies"] = model
        models["prot"][i]["bins"] = edges
        model, edges = get_null_roll(sequences, i)
        models["roll"][i] = {}
        models["roll"][i]["frequencies"] = model
        models["roll"][i]["bins"] = edges
        model, edges = get_null_helt(sequences, i)
        models["helt"][i] = {}
        models["helt"][i]["frequencies"] = model
        models["helt"][i]["bins"] = edges
        print(models)
        print("Finsihed calculating null for shapes of length {}".format(i))
        with open("models/models", "wb") as outfile:
            pickle.dump(models, outfile)
