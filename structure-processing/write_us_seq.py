#/usr/bin/env python2
import sys
import numpy as np
import pickle
import matplotlib.pyplot as plt
import os
from collections import Counter
import complicated_heatmap
#import pandas as pd

data_path = '/Users/anatale/Documents/school/UCSF/Kortemme_lab/data/second_run_scores'
# with open(os.path.join(data_path, 'backrub_scores.pkl'), 'rb') as backrub_scores_file:
#     backrub_scores = pickle.load(backrub_scores_file)

fixbb_scores = []
for state in [1,2,3,4]:
    with open(os.path.join(data_path, 'state_%s_data_v2.pkl' % state), 'rb') as data_file:
        data_array, data_labels = pickle.load(data_file)
    fixbb_scores.append((data_array,data_labels))

yeast_gsp1_seq = '\
MSAPAANGEVPTFKLVLVGDGGTGKTTFVKRHLTGEFEKKYIATIGVEVHPLSFYTNFGE\
IKFDVWDTAGQEKFGGLRDGYYINAQCAIIMFDVTSRITYKNVPNWHRDLVRVCENIPIV\
LCGNKVDVKERKVKAKTITFHRKKNLQYYDISAKSNYNFEKPFLWLARKLAGNPQLEFVA\
SPALAPPEVQVDEQLMQQYQQEMEQATALPLPDEDDADL'

# in array:
# 1st index is substate id
# 2nd index is microstate id
# 3nd index is mutation
# 4th index is residue number

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

numbertoamino = {}
for elem in aminotonumber:
    numbertoamino[aminotonumber[elem]] = elem

mutation_labels = []
for value in sorted(aminotonumber.values()):
    for key in aminotonumber:
        if key != '*':
            if value == aminotonumber[key]:
                mutation_labels.append(key)
mutation_labels.reverse()

sequence_labels = ['K101', 'N102', 'V103', 'N105', 'W106', 'H107', 'R108', 'D109', 'C122','G123','N124','K125','V126','D127']

state_names = ['apo +GEF','GTP +GAP','GDP','GTP +importin']

wt_indicators = [aminotonumber[s[0]]-1 for s in sequence_labels]

with open('fixbb_seqs.txt', 'w') as outfile:
    outfile.write('positions\n')
    for position in sequence_labels:
        outfile.write('%s ' % position[1:])
    outfile.write('\n')
    outfile.write('native res\n')
    for position in sequence_labels:
        outfile.write('%s' % position[0])
    outfile.write('\n')
    for state_idx, name in enumerate(state_names):
        state_labels = fixbb_scores[state_idx][1]
        state_array = fixbb_scores[state_idx][0]
        outfile.write('%s\n' % name)
        print name
        for sub_label in sorted(state_labels):
            sub_idx = state_labels[sub_label][0]
            outfile.write('%s\n' % sub_label)
            print sub_label
            for us_name in sorted(state_labels[sub_label][1]):
                us_idx = state_labels[sub_label][1][us_name]
                if len(str(us_name)) == 1:
                    outfile.write('%s   ' % str(us_name))
                elif len(str(us_name)) == 2:
                    outfile.write('%s  ' % str(us_name))
                elif len(str(us_name)) == 3:
                    outfile.write('%s ' % str(us_name))
                print us_name
                for pos_idx in range(np.size(state_array,3)):
                    #print np.min(state_array[sub_idx,us_idx,:,pos_idx])
                    print numbertoamino[state_array[sub_idx,us_idx,:,pos_idx].argmin() + 1]
                    outfile.write('%s' % numbertoamino[state_array[sub_idx,us_idx,:,pos_idx].argmin() + 1])
                outfile.write('\n')
