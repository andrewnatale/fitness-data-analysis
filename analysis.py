#!/usr/bin/env python2
#
# analyze gsp1 fitness data
"""
what this program does:

    - loads wt sequence data into a dictionary using fasta_aa module

    - loads sequence enrichment data from a csv file with
    the format "position,aa,v1,v2,...,vn" where each v is an experimental
    repeat

    - turns this data into a dictionary with the format
    (position, aa, #) as keys and the mean of v1 and v2 as values
    '#' is a code for each aa used to order them by chemical properties

    - normalizes the fitness data where zero is the average stop codon
    enrichment and 1 is the average wt fitness. values below zero are set
    to zero

what this program needs to do:
    - make a better heatmap using normalized data
    - get data into a format to enter into weblogo
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
first_resi = sorted([i[0] for i in fitness])[0]

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
            if elem in self.bad_data:
                pass
            else:
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
        self.hmap_array = np.zeros((21,seq_length))
        # fill it with data
        for elem in self.fitness_norm:
                hmap_array[int(elem[2]), int(elem[0]) - first_resi] = \
                    fitness_normalized[data]


# for each experimental set, build an object and put it in a list
fitness_sets = []
for i in range(expreps):
    fitness_set = {}
    for elem in fitness:
        fitness_set[elem] = fitness[elem][i]
    fitness_sets.append(ExpReps(fitness_set))

for rep in fitness_sets:
    rep.normalize()
    for elem in sorted(rep.fitness_norm.values()):
        print elem

# compose heatmap axis labels
row_labels = []
for n in range(101, 141, 1):
    row_labels.append(('%s%d') % (native_seq[n], n))

column_labels = []
for value in sorted(aminotonumber.values()):
    for key in aminotonumber:
        if value == aminotonumber[key]:
            column_labels.append(key)

#plot_hmap(hmap_array, row_labels, column_labels)

"""
# extract stop codon data points into a new dict
stop_codon_fitness = {}
for elem in fitness:
    if elem[1] == '*':
        stop_codon_fitness[elem] = fitness[elem]

# get avg enrichment of stop codon mutants
stop_avg = sum(stop_codon_fitness.values()) / len(stop_codon_fitness.values())



# get the average of wt positions in the shifted data
sum_wt = 0
count_wt = 0
for elem in fitness_shifted:
    for wt_pos in native_seq:
        if int(elem[0]) == wt_pos and elem[1] == native_seq[wt_pos]:
            sum_wt += fitness_shifted[elem]
            count_wt += 1
avg_shifted_wt = sum_wt / count_wt
#print stop_avg
#print avg_shifted_wt
# now divide every shifted point by the wt average
fitness_normalized = {}
for elem in fitness_shifted:
    fitness_normalized[elem] = fitness_shifted[elem] / avg_shifted_wt

# set negative values to zero
for elem in fitness_normalized:
    if fitness_normalized[elem] < 0:
        fitness_normalized[elem] = 0
#    if fitness_normalized[elem] > 100:
#        fitness_normalized[elem] = 0

# make an empty array the size of the heatmap
hmap_array = np.zeros((21,40))

# fill it with data
for data in fitness_normalized:
        position = int(data[0])
        mutation = int(data[2])
        hmap_array[mutation, position - 101] = fitness_normalized[data]

# compose heatmap axis labels
row_labels = []
for n in range(101, 141, 1):
    row_labels.append(('%s%d') % (native_seq[n], n))

column_labels = []
for value in sorted(aminotonumber.values()):
    for key in aminotonumber:
        if value == aminotonumber[key]:
            column_labels.append(key)

#plot_hmap(hmap_array, row_labels, column_labels)
"""
