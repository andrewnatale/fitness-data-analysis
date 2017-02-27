import sys
import pandas as pd
import os
import numpy as np
import pickle
import csv

scorefile_dir = '/Users/anatale/school/UCSF/Kortemme_lab/data/combo_scores/combined'
#output_dir = os.getcwd()

seq_index = pd.read_table('/Users/anatale/school/UCSF/Kortemme_lab/code/multi-state-design/structure-processing/aux_files/Gsp1_index_all.txt')

yeast_gsp1_seq = '\
MSAPAANGEVPTFKLVLVGDGGTGKTTFVKRHLTGEFEKKYIATIGVEVHPLSFYTNFGE\
IKFDVWDTAGQEKFGGLRDGYYINAQCAIIMFDVTSRITYKNVPNWHRDLVRVCENIPIV\
LCGNKVDVKERKVKAKTITFHRKKNLQYYDISAKSNYNFEKPFLWLARKLAGNPQLEFVA\
SPALAPPEVQVDEQLMQQYQQEMEQATALPLPDEDDADL'

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


# macrostates - dict used for organizing other data structures
macrostates = {
1 : [
('1I2M_TP0_backrub', [1,2,3,4,5,6,7,8])
],
2 : [
('1A2K_TP0_backrub', [1,2,3,4,5,6,7,8,9,10])
],
3 : [
('1K5D_ALL_backrub', [1,2,3,4,5])
],
4 : [
('3GJ0_TP0_backrub', [1,2,3,4,5,6])
],
5 : [
('1QBK_TP0_backrub', [1,2]),
('2BKU_TP0_backrub', [1,2]),
('3W3Z_TP0_backrub', [1,2]),
('4OL0_TP0_backrub', [1,])
],
6 : [
('1I2M_NON_backrub', [1,]),
('1K5D_NON_backrub', [1,]),
('1QBK_NON_backrub', [1,]),
('3GJ0_TP0_backrub', [7,]),
('1A2K_NON_backrub', [1,])
]
}

target_pos = [101,102,103,104,105,106,107,108,109,122,123,124,125,126,127]
seq_to_arr_idx = {}
for idx,pos in enumerate(target_pos):
    seq_to_arr_idx[pos] = idx

# turn the index file into a dictionary for much faster lookups
seq_dict = {}
for index, row in seq_index.iterrows():
    seq_dict[(row['pdb'][0:4], row['pdb_res_num'])] = row['yeast_seq_num']

# given a structure_id and a res number in pdb numbering, return the yeast seq number and wt aa id
def lookup_res(structure_id, pdb_num):
    try:
        yeast_num = seq_dict[(structure_id.lower(), pdb_num)]
        wt_aa_id = yeast_gsp1_seq[yeast_num-1]
        return yeast_num, wt_aa_id
    except KeyError:
        print('structure id not found in index file')
        raise

def parse_state(state_id):
    # unpack score files and do inital counting
    scores_for_state = {}
    for substate in macrostates[state_id]:
        for i in substate[1]:
            substate_name = '%s_%s' % (substate[0], i)
            if os.path.isfile(os.path.join(scorefile_dir, '%s+.sc' % substate_name)):
                print('processing ', substate_name)
            else:
                print('score file for %s not found, exiting' % substate_name)
                exit(1)
            with open(os.path.join(scorefile_dir, '%s+.sc' % substate_name), 'r') as data_file:
                all_lines = data_file.readlines()
            data_lines = [i.strip('\n').split() for i in all_lines[2:]]
            # figure out which microstates are used
            ensemble_ids = set()
            for line in data_lines:
                iden = line[-1].split('_')
                ensemble_ids.add(int(iden[2][7:]))
            # store the extracted data in a dictionary
            scores_for_state[substate_name] = (sorted(ensemble_ids), data_lines)
    # init the array, construct labels
    state_size = 0
    state_labels = {}
    for substate in scores_for_state:
        #state_size += len(scores_for_state[substate][0])
        for microstate_id in scores_for_state[substate][0]:
            structure_name = '_'.join(substate.split('_')[0:-1])
            state_labels[(structure_name, microstate_id)] = state_size
            state_size += 1
    # state array dimensions
    # 1) # of microstates
    # 2) # of amino acids (fixed at 20)
    # 3) # of positions considered
    state_array = np.zeros((state_size, 20, len(target_pos)), dtype=np.float64)
    # fill array
    for substate in scores_for_state:
        for line in scores_for_state[substate][1]:
            #print(line)
            score = float(line[1])
            position = int(line[-1].split('_')[-1][:-1])
            mutation = line[-1][-1]
            structure_name = line[-1][:16]
            microstate_id = int(line[-1].split('_')[2][7:])
            us_idx = state_labels[(structure_name, microstate_id)]
            mut_idx = aminotonumber[mutation] - 1
            pos_idx = seq_to_arr_idx[lookup_res(structure_name[0:4], position)[0]]
            #print(us_idx,mut_idx,pos_idx)
            #print(structure_name, microstate_id, position, mutation, score)
            state_array[us_idx,mut_idx,pos_idx] = score
    # print(state_array)
    # print(state_array.shape)
    return (state_array, state_labels)

for state in macrostates:
    stuff = parse_state(state)
    pickle.dump(stuff, open('state_%s_data_v5.pkl' % str(state), 'wb'))

# # parse and organize all the score files
# score_data = {}
# for key in macrostates:
#     for substate in macrostates[key]:
#         for i in substate[1]:
#             substate_name = '%s_%s' % (substate[0], i)
#             if os.path.isfile(os.path.join(scorefile_dir, '%s+.sc' % substate_name)):
#                 # store a tuple for each substate, where the first value is the macrostate id,
#                 # the second is a list of the backrub ensemble ids that went into the substate,
#                 # and third is the dataframe containing scores for the substate
#                 with open(os.path.join(scorefile_dir, '%s+.sc' % substate_name), 'r') as data_file:
#                     data_lines = data_file.readlines()
#                 data_lines = [i.strip('\n').split() for i in data_lines[2:]]
#                 ensemble_ids = set()
#                 for line in data_lines:
#                     iden = line[-1].split('_')
#                     ensemble_ids.add(int(iden[2][7:]))
#                 score_data[substate_name] = (key, sorted(ensemble_ids), data_lines)

# for key in sorted(score_data):
#     print key, score_data[key][0], score_data[key][1], len(score_data[key][2])
# print len(score_data)

# # organize all the scores for a specified macrostate into an array for processing
# def build_mstate_array(macrostate_id):
#     # build dict for indexing data array
#     datalabels = {}
#     datalabel_names = []
#     for substate in score_data:
#         if score_data[substate][0] == macrostate_id:
#             datalabel_names.append(substate)
#     for idx,substate in enumerate(sorted(datalabel_names)):
#         datalabels[substate] = (idx, {})
#         for jdx,ustate_id in enumerate(score_data[substate][1]):
#             datalabels[substate][1][ustate_id] = jdx
#
#     # for datalabel in sorted(datalabels):
#     #     print datalabel, datalabels[datalabel][0]
#     #     for us_id in sorted(datalabels[datalabel][1]):
#     #         print us_id, datalabels[datalabel][1][us_id]
#
#     # build an empty array to hold data
#     # 1st index is substate id
#     # 2nd index is microstate id
#     # 3nd index is mutation
#     # 4th index is residue number
#     datarray = np.zeros((len(datalabels), 10, 20, len(target_pos)))
#     # fill the array up, iterating over each substate
#     for ss_id in datalabels:
#         # get corresponding data
#         print('processing %s scores' % ss_id)
#         for line in score_data[ss_id][2]:
#             iden = line[-1].split('_')
#             rel_seq_idx, nataa = lookup_res('%s_%s' % (iden[0], iden[1]), int(iden[4][:-1]))
#             #print rel_seq_idx, nataa
#             print('%s_%s' % (iden[0], iden[1]), iden[4][:-1], iden[4][-1], line[1], rel_seq_idx, nataa)
#             # print iden
#             # insert score into the array
#             # print datalabels[ss_id][0]
#             # print datalabels[ss_id][1][int(iden[2][7:])]
#             # print aminotonumber[iden[4][-1]]-1
#             # print seq_to_arr_idx[int(rel_seq_idx)]
#             datarray[datalabels[ss_id][0],
#                      datalabels[ss_id][1][int(iden[2][7:])],
#                      aminotonumber[iden[4][-1]]-1,
#                      seq_to_arr_idx[int(rel_seq_idx)]] = float(line[1])
#     return (datarray, datalabels)
#
# for state in macrostates:
#     #print build_mstate_array(state)
#     state_data_array = build_mstate_array(state)
#     pickle.dump(state_data_array, open('state_%s_data_v4.pkl' % str(state), 'wb'))
