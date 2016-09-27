#!/usr/bin/env python2
#
# analyze gsp1 fitness data
#
# import this module to process fitness data from a csv file with the format
# "position,aa,v1,v2,...,vn" where each v is an experimental repeat and
# return that data with some useful labels for visualization

# imports
import numpy as np
from collections import Counter
from fasta_aa import load_fasta_aa

# main function, import this into other scripts
def process_fitness(fitness_data_file, native_seq_file):

    # mapping of aa types to numbers that orders them from
    # most hydrophobic to most polar
    aminotonumber = {'*': 0, 'A': 9, 'C': 8, 'E': 20, 'D': 19,
                     'G': 10, 'F': 2, 'I': 5, 'H': 16, 'K': 18, 'M': 6,
                     'L': 4, 'N': 14, 'Q': 15, 'P': 11, 'S': 12, 'R': 17,
                     'T': 13, 'W': 1, 'V': 7, 'Y': 3}

    # load a reference sequence - this will be used to normalize the data
    # as well as to provide labels
    native_seq = load_fasta_aa(native_seq_file)

    # open and cleanup enrichment data
    data = open(fitness_data_file, 'r')
    datalines = data.readlines()
    data.close()
    datalines = [i.strip('\r\n') for i in datalines]
    # remove column labels
    fitness_header = datalines.pop(0)

    # count the number of replicates in the input file
    expreps = len(datalines[0].split(',')) - 2

    # read data into a dictionary, replacing the aa leter code with a number
    # which will be used to arrange them (see 'aminotonumber' dictionary)
    fitness_dict = {}
    for line in datalines:
        elem = line.split(',')
        fitness_dict[(elem[0], aminotonumber[elem[1]])] = \
          [float(i) for i in elem[2:]]

    # get length of sequence and first resi position
    seq_length = len(Counter([i[0] for i in fitness_dict]))
    first_resi = int(sorted([i[0] for i in fitness_dict])[0])

    # create an empty 3d array of right size for all the data
    fitness_array = np.zeros((expreps,21,seq_length))

    # for each replicate, fill the array with data and process it
    for i in range(expreps):

        # fill the array using the data dictionary
        for elem in fitness_dict:
            fitness_array[i, int(elem[1]), int(elem[0]) - first_resi] \
              = fitness_dict[elem][i]

        # look for bad data points within this replicate
        bad_data_exprep = []
        for index, j in np.ndenumerate(fitness_array[i,:,:]):
            if j == 999.9:
                # remember that data point, then set it to zero for convenience
                bad_data_exprep.append(index)
                fitness_array[i,index[0],index[1]] = 0.0

        # get average stop codon values for this replicate
        avg_stop = np.sum(fitness_array[i,aminotonumber['*'],:]) \
          / fitness_array[i,aminotonumber['*'],:].__len__()
        # shift everything by the average stop codon value
        fitness_array[i,:,:] = fitness_array[i,:,:] - avg_stop

        # get the average wt residue fitness in the shifted data
        avg_wt = 0
        for j in range(first_resi, first_resi + seq_length, 1):
            if (aminotonumber[native_seq[j]],j - first_resi) \
              not in bad_data_exprep:
                avg_wt += \
                  fitness_array[i,aminotonumber[native_seq[j]],j - first_resi]
        avg_wt = avg_wt / (seq_length - len(bad_data_exprep))
        # divide everything by the average wt residue fitness
        fitness_array[i,:,:] = fitness_array[i,:,:] / avg_wt

        # set any values below zero to zero
        for index, j in np.ndenumerate(fitness_array[i,:,:]):
            if j < 0:
                fitness_array[i,index[0],index[1]] = 0.0

    # compose axis labels
    column_labels = []
    for n in range(first_resi, first_resi + seq_length, 1):
        column_labels.append(('%s%d') % (native_seq[n], n))

    row_labels = []
    for value in sorted(aminotonumber.values()):
        for key in aminotonumber:
            if value == aminotonumber[key]:
                row_labels.append(key)

    return (np.mean(fitness_array, axis=0), column_labels, row_labels)
