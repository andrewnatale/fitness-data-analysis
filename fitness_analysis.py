#!/usr/bin/env python2

# process fitness data from a csv file with the format
# "position,aa,v1,v2,...,vn" where each v is an experimental repeat and
# return that data with some useful labels for visualization.
# also requires a native sequence file in fasta format for normalization

# normalization types:
# none -> None
# zero(max stop) only -> 'stop'
# per position probability -> 'prob'

import numpy as np
from collections import Counter
from Bio import SeqIO
import matplotlib.pyplot as plt
import math
import complicated_heatmap

normalization = 'stop'

# hardcoded data filenames
native_seq_file = '../../notebook/logo/P32835.fasta'
fitness_data_file = '../../notebook/logo/gsp_log_data_101-140.csv'

# mapping of aa types to numbers that loosely
# orders them from hydrophobic to polar with 0 as 'stop'
aminotonumber = {'*': 0, 'W': 1,
                 'F': 2, 'Y': 3,
                 'L': 4, 'I': 5,
                 'M': 6, 'V': 7,
                 'C': 8, 'A': 9,
                 'G': 10, 'P': 11,
                 'S': 12, 'T': 13,
                 'N': 14, 'Q': 15,
                 'H': 16, 'R': 17,
                 'K': 18, 'D': 19,
                 'E': 20}

# normalization types:
# None: don't do any normalization, results in a mix of pos and neg values
def norm_none(fitness_array):
    return np.mean(fitness_array, axis=0)

# 'stop': use the highest stop codon fitness as a cutoff, set it to
# zero and shift all the data, then discard values below zero
def norm_stop(fitness_array):
    mean_fitness_array = np.mean(fitness_array, axis=0)
    # get max stop codon values
    max_stop = np.amax(mean_fitness_array[aminotonumber['*'],:])
    # shift everything by the max stop codon value
    mean_fitness_array = mean_fitness_array - max_stop
    # set bad data points and values below zero to zero
    for index, j in np.ndenumerate(mean_fitness_array):
        if j < 0:
            mean_fitness_array[index] = 0.0
    # get the average wt residue fitness in the shifted data
    wt_fitness = np.empty(0)
    for j in range(first_resi, first_resi + seq_length, 1):
        wt_fitness = np.append(wt_fitness, mean_fitness_array[aminotonumber[native_seq[0].seq[j-1]],j - first_resi])
    wt_avg = np.mean(wt_fitness[~np.isnan(wt_fitness)])
    # divide everything by the average wt residue fitness
    return mean_fitness_array / wt_avg

# # 'prob': same as stop, but then normalize each position column
# # so that it sums to 1
# elif normalization == 'prob':
#     for i in range(expreps):
#         # get max stop codon values for this replicate
#         max_stop = np.amax(fitness_array[i,aminotonumber['*'],:])
#         # shift everything by the max stop codon value
#         fitness_array[i,:,:] = fitness_array[i,:,:] - max_stop
#     # set any values below zero to zero
#     for index, j in np.ndenumerate(fitness_array):
#         if j < 0:
#             fitness_array[index] = 0.0
#         if index in bad_data:
#             fitness_array[index] = np.nan
#     # average all the expreps
#     mean_fitness_array = np.mean(fitness_array, axis=0)
#     # normalize each position independently
#     for i in range(seq_length):
#         mean_fitness_array[:,i] = mean_fitness_array[:,i] \
#           / np.sum(mean_fitness_array[:,i])


# load a reference sequence - this will be used to normalize the data
# as well as to provide labels. this should be a fasta file and only
# the first record in the file will be used
native_seq = list(SeqIO.parse(native_seq_file, 'fasta'))

# open fitness data
data = open(fitness_data_file, 'r')
datalines = data.readlines()
data.close()
datalines = [i.strip('\r\n') for i in datalines]
# remove column labels (but save it just in case)
fitness_header = datalines.pop(0)

# count the number of replicates in the input file
expreps = len(datalines[0].split(',')) - 2

# read data into a dictionary, replacing the aa letter code with a number
# which will be used to arrange them (see 'aminotonumber' dictionary)
fitness_dict = {}
for line in datalines:
    elem = line.split(',')
    fitness_dict[(elem[0], aminotonumber[elem[1]])] = [float(i) for i in elem[2:]]

# get length of sequence and first resi position
seq_length = len(Counter([i[0] for i in fitness_dict]))
first_resi = int(sorted([i[0] for i in fitness_dict])[0])

# create an empty 3d array of the correct size for all the data
fitness_array = np.zeros((expreps,21,seq_length))

# for each replicate, fill the array with data and process it.
for i in range(expreps):
    # fill the array using the data dictionary
    for elem in fitness_dict:
        fitness_array[i, int(elem[1]), int(elem[0]) - first_resi] = fitness_dict[elem][i]
# look for bad data points
print 'Bad data points:'
for index, j in np.ndenumerate(fitness_array):
    if j == 999.9:
        fitness_array[index] = np.nan
        print index

if normalization == None:
    processed_array = norm_none(fitness_array)
elif normalization == 'stop':
    processed_array = norm_stop(fitness_array)
elif normalization == 'prob':
    pass

# compose data labels
sequence_labels = []
for n in range(first_resi, first_resi + seq_length, 1):
    sequence_labels.append(('%s%d') % (native_seq[0].seq[n-1], n))

mutation_labels = []
for value in sorted(aminotonumber.values()):
    for key in aminotonumber:
        if value == aminotonumber[key]:
            mutation_labels.append(key)
mutation_labels.reverse()

complicated_heatmap.plot_map(processed_array[1:,21:28], mutation_labels[:-1], sequence_labels[21:28], 1, True)
