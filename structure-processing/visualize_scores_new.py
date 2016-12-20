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

# adjust all scores to be relative to the wt aa score at that position
def fixbb_diff(state):
    state_idx = state - 1
    state_labels = fixbb_scores[state_idx][1]
    state_array = fixbb_scores[state_idx][0]
    for sub_label in sorted(state_labels):
        #structure_name = substate_label[:-10]
        sub_idx = state_labels[sub_label][0]
        #print substate_label, substate_index, structure_id
        for us_name in sorted(state_labels[sub_label][1]):
            us_idx = state_labels[sub_label][1][us_name]
            #print us_id, us_index
            #print fixbb_scores[state-1][0][substate_index,us_index,:,:]
            for pos_idx, mut_idx in enumerate(wt_indicators):
                state_array[sub_idx,us_idx,:,pos_idx] -= state_array[sub_idx,us_idx,mut_idx,pos_idx]
            #print fixbb_scores[state-1][0][substate_index,us_index,:,:]
    return state_array

scale_factor = 2.0
state_listA = []
state_listB = []
for state in [1,2,3,4]:
    state_idx = state - 1
    state_array = fixbb_diff(state)
    avg_state_array = np.mean(state_array, (0,1))
    exp_array = np.exp((-1.0 * avg_state_array) / scale_factor)
    state_listA.append(exp_array[:,8:])
    state_listB.append(exp_array[:,0:8])
    #complicated_heatmap.single_map(exp_array[:,8:], mutation_labels, sequence_labels, use_mpoint=True, mpoint=1.0, set_limits=(2.,0.0), show_cbar=True)
complicated_heatmap.new_multi_map(state_listA, mutation_labels, sequence_labels[8:], desc_list=state_names, use_mpoint=True, mpoint=1.0, set_limits=(2.0,0.0), show_cbar=True)
complicated_heatmap.new_multi_map(state_listB, mutation_labels, sequence_labels[0:8], desc_list=state_names, use_mpoint=True, mpoint=1.0, set_limits=(2.0,0.0), show_cbar=True)

# transform scores to differences
# wt_indicators = [aminotonumber[s[0]]-1 for s in sequence_labels]
# for state in [1,2,3,4]:
#     calc_fixbb_diff(state, wt_indicators)

# for state in [1,2,3,4,5]:
#     calc_difference(state)

# stdev = False
#
# if stdev:
#     new_state_list = []
#     for state in [1,2,3,4]:
#        state_std = np.std(fixbb_scores[state-1][0], (0,1))
#        new_state_list.append(state_std[:,8:])
#     complicated_heatmap.new_multi_map(new_state_list, mutation_labels, sequence_labels[8:], desc_list=state_names, set_limits=(200.,0.), show_cbar=True)
# else:
#     new_state_list = []
#     for state in [1,2,3,4]:
#        state_avg = np.mean(fixbb_scores[state-1][0], (0,1))
#        new_state_list.append(state_avg[:,0:8])
#     complicated_heatmap.new_multi_map(new_state_list, mutation_labels, sequence_labels[0:8], desc_list=state_names, use_mpoint=True, set_limits=(15.,-5.), show_cbar=True)

# a = 0
# state_a_avg = np.mean(fixbb_scores[a][0], (0,1))
# b = 3
# state_b_avg = np.mean(fixbb_scores[b][0], (0,1))
#
# delta_ab = state_a_avg - state_b_avg
#
# toggle = False
#
# if toggle:
#     complicated_heatmap.single_map(delta_ab[:,0:8], mutation_labels, sequence_labels[0:8], use_mpoint=True, set_limits=(50.,-50.), show_cbar=True)
# else:
#     complicated_heatmap.single_map(delta_ab[:,8:], mutation_labels, sequence_labels[8:], use_mpoint=True, set_limits=(50.,-50.),show_cbar=True)
