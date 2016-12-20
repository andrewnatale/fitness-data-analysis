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

#print data_labels

#print backrub_scores
def fixbb_hist(scale=None):
    fig = plt.figure()
    a = fig.add_subplot(2,3,1)
    plt.hist(fixbb_scores[0][0].flatten(), bins=150, range=scale)
    a.set_title('macrostate 1')
    a = fig.add_subplot(2,3,2)
    plt.hist(fixbb_scores[1][0].flatten(), bins=150, range=scale)
    a.set_title('macrostate 2')
    a = fig.add_subplot(2,3,3)
    plt.hist(fixbb_scores[2][0].flatten(), bins=150, range=scale)
    a.set_title('macrostate 3')
    a = fig.add_subplot(2,3,4)
    plt.hist(fixbb_scores[3][0].flatten(), bins=150, range=scale)
    a.set_title('macrostate 4')
    # a = fig.add_subplot(2,3,5)
    # plt.hist(fixbb_scores[4][0].flatten(), bins=150, range=scale)
    # a.set_title('macrostate 5')
    plt.show()

# raw_scores_hist()

def backrub_hist():
    structure_ids = sorted(Counter([i[0] for i in backrub_scores]).keys())
    fig = plt.figure()
    plot_id = 1
    for structure_id in structure_ids:
        structure_backrub_scores = [backrub_scores[key] for key in backrub_scores if key[0] == structure_id]
        a = fig.add_subplot(3,6,plot_id)
        plt.hist(structure_backrub_scores, bins=50)
        a.set_title(structure_id)
        plot_id += 1
    plt.show()

# backrub_hist()

# for each fixbb score, subtract the score of the backrub structure from which it originated
# then plot the distribution of the differences
# this transforms the array in place
def calc_difference(state):
    for substate_label in sorted(fixbb_scores[state-1][1]):
        structure_id = substate_label[:-10]
        substate_index = fixbb_scores[state-1][1][substate_label][0]
        #print substate_label, substate_index, structure_id
        for us_id in sorted(fixbb_scores[state-1][1][substate_label][1]):
            us_index = fixbb_scores[state-1][1][substate_label][1][us_id]
            #print us_id, us_index
            backrub_score = backrub_scores[(structure_id, str(us_id))]
            #print fixbb_scores[state-1][0][substate_index,us_index,:,:]
            #print backrub_score
            fixbb_scores[state-1][0][substate_index,us_index,:,:] -= backrub_score
            #print fixbb_scores[state-1][0][substate_index,us_index,:,:]

# similar to above, except use the fixbb structure with the wt res as the
# baseline score to subtract
def calc_fixbb_diff(state, wt_indicators):
    for substate_label in sorted(fixbb_scores[state-1][1]):
        structure_id = substate_label[:-10]
        substate_index = fixbb_scores[state-1][1][substate_label][0]
        #print substate_label, substate_index, structure_id
        for us_id in sorted(fixbb_scores[state-1][1][substate_label][1]):
            us_index = fixbb_scores[state-1][1][substate_label][1][us_id]
            #print us_id, us_index
            #print fixbb_scores[state-1][0][substate_index,us_index,:,:]
            for idx, elem in enumerate(wt_indicators):
                wt_fixbb_score = fixbb_scores[state-1][0][substate_index,us_index,elem,idx]
                fixbb_scores[state-1][0][substate_index,us_index,:,idx] -= wt_fixbb_score
            #print fixbb_scores[state-1][0][substate_index,us_index,:,:]

#fixbb_hist()

#fixbb_hist(scale=[-100, 100])

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

mutation_labels = []
for value in sorted(aminotonumber.values()):
    for key in aminotonumber:
        if key != '*':
            if value == aminotonumber[key]:
                mutation_labels.append(key)
mutation_labels.reverse()
#print mutation_labels

sequence_labels = ['C122','G123','N124','K125','V126','D127']

state_names = ['apo +GEF','GTP +GAP','GDP','GTP +importin','various']

# transform scores to differences
wt_indicators = [aminotonumber[s[0]]-1 for s in sequence_labels]
for state in [1,2,3,4,5]:
    calc_fixbb_diff(state, wt_indicators)

# for state in [1,2,3,4,5]:
#     calc_difference(state)

new_state_list = []
for state in [1,2,3,4,5]:
    state_avg = np.mean(fixbb_scores[state-1][0], (0,1))
    new_state_list.append(state_avg)
complicated_heatmap.multi_map(new_state_list, state_names, mutation_labels, sequence_labels, 0, set_limits=(20.,-20.), show_cbar=True)
