#!/usr/bin/env python2
#
# analyze gsp1 fitness data
"""
what this program does:

    - loads baseline sequence data into a dictionary using fasta_aa module

    - takes mutation enrichment data from a csv file with the format
    "position,aa,v1,v2,...,vn" where each v is an experimental repeat

    - normalizes the fitness data where zero is the average stop codon
    enrichment and 1 is the average baseline fitness; values below zero
    are set to zero

    - plots a heatmap of the normalized data

what this program needs to do:

    - use a three dimensional array or dataframe instead of a bunch
    of dictionaries

    - get data into a format to enter into weblogo to generate a sequence
    profile logo
"""

import numpy as np
from collections import Counter
from fasta_aa import load_fasta_aa
from heatmap import plot_hmap

# hardcoded data filenames
native_seq_file = 'P32835.fasta'
fitness_data_file = 'gsp_log_data_101-140.csv'

# mapping of aa types to numbers that orders them from
# most hydrophobic to most polar
aminotonumber = {'*': 0, 'A': 9, 'C': 8, 'E': 20, 'D': 19,
                 'G': 10, 'F': 2, 'I': 5, 'H': 16, 'K': 18, 'M': 6,
                 'L': 4, 'N': 14, 'Q': 15, 'P': 11, 'S': 12, 'R': 17,
                 'T': 13, 'W': 1, 'V': 7, 'Y': 3}

native_seq = load_fasta_aa(native_seq_file)

# open and cleanup enrichment data
data = open(fitness_data_file, 'r')
datalines = data.readlines()
data.close()
datalines = [i.strip('\r\n') for i in datalines]
# remove column labels
fitness_header = datalines.pop(0)

expreps = len(datalines[0].split(',')) - 2

# read data into a dictionary
fitness = {}
for line in datalines:
    elem = line.split(',')
    fitness[(elem[0], elem[1], aminotonumber[elem[1]])] = \
      [float(i) for i in elem[2:]]

# length of sequence and first resi position
seq_length = len(Counter([i[0] for i in fitness]))
first_resi = int(sorted([i[0] for i in fitness])[0])

class ExpReps(object):
    def __init__(self, exprep):
        self.fitness = exprep
        # there are some wierd entries where the fitness value is 999.9
        # set these to 0 and remember their names for later
        self.bad_data = []
        for elem in self.fitness:
            if self.fitness[elem] == 999.9:
                self.bad_data.append(elem)
                self.fitness[elem] = 0.0

    def normalize(self):
        # get avg enrichment of stop codon mutants
        self.avg_stop = 0
        self.stop_count = 0
        for elem in self.fitness:
            if elem[1] == '*':
                self.avg_stop += float(self.fitness[elem])
                self.stop_count += 1
        self.avg_stop = self.avg_stop / float(self.stop_count)
        # start normalizing by shifting eveything up by the magnitude
        # of the avg STOP codon fitness (want this to be our zero)
        self.fitness_norm = {}
        for elem in self.fitness:
            self.fitness_norm[elem] = self.fitness[elem] - self.avg_stop
        # get the average of wt positions in the shifted data
        self.sum_wt = 0
        self.count_wt = 0
        for elem in self.fitness_norm:
            if elem not in self.bad_data:
                for wt_pos in native_seq:
                    if int(elem[0]) == wt_pos \
                      and elem[1] == native_seq[wt_pos]:
                        self.sum_wt += self.fitness_norm[elem]
                        self.count_wt += 1
        self.avg_shifted_wt = self.sum_wt / self.count_wt
        # now divide every shifted point by the wt average, then clip
        # negative values to 0
        for elem in self.fitness_norm:
            self.fitness_norm[elem] = self.fitness_norm[elem] \
              / self.avg_shifted_wt
            if self.fitness_norm[elem] < 0.0:
                self.fitness_norm[elem] = 0.0

    def build_array(self):
        # make an empty array the size of the heatmap
        self.hmap_array = np.zeros((20,seq_length))
        # fill it with data
        for elem in self.fitness_norm:
            if elem[1] != '*':
                self.hmap_array[int(elem[2] - 1), int(elem[0]) - first_resi] \
                  = self.fitness_norm[elem]


# for each experimental set, build an object and put it in a list
fitness_reps = []
for i in range(expreps):
    fitness_rep = {}
    for elem in fitness:
        fitness_rep[elem] = fitness[elem][i]
    fitness_reps.append(ExpReps(fitness_rep))

for rep in fitness_reps:
    rep.normalize()
    rep.build_array()

# average the experimental sets
hmap_norm_mean = np.zeros((20,seq_length))
for i in range(expreps):
    hmap_norm_mean += fitness_reps[i].hmap_array
hmap_norm_mean = hmap_norm_mean / float(expreps)

# compose heatmap axis labels
row_labels = []
for n in range(first_resi, first_resi + seq_length, 1):
    row_labels.append(('%s%d') % (native_seq[n], n))

column_labels = []
for value in sorted(aminotonumber.values()):
    for key in aminotonumber:
        if value == aminotonumber[key]:
            column_labels.append(key)

plot_hmap(hmap_norm_mean, row_labels, column_labels)
