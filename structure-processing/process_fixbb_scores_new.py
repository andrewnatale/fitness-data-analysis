import sys
import pandas as pd
import os
import numpy as np
import pickle
import csv

scorefile_dir = '/Users/anatale/Documents/school/UCSF/Kortemme_lab/data/second_run_scores'

seq_index = pd.read_table('/Users/anatale/Documents/school/UCSF/Kortemme_lab/code/structure-processing/aux_files/Gsp1_index_all.txt')

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
# GEF-complexed, no nucleotide - one structure: 1I2M
1 : [
('1I2M_TP0_backrub', [1,2,3,4,5,6])
],
# GAP-complexed, GTP-analogue nucleotide - one structure: 1K5D
2 : [
('1K5D_TP0_backrub', [1,2,3,4,5,6])
],
# uncomplexed, GDP nucleotide - three structures: 3GJ0, 5CIQ, 5CIT
3 : [
('3GJ0_TP0_backrub', [1,2]),
('5CIQ_TP0_backrub', [1,2]),
('5CIT_TP0_backrub', [1,2])
],
# karyopherin complexed, GTP nucleotide - six structures: 1QBK, 1WA5_TP2, 2BKU, 3ICQ, 3M1I_TP2, 3W3Z
4 : [
('1QBK_TP0_backrub', [1,]),
('1WA5_TP2_backrub', [1,]),
('2BKU_TP0_backrub', [1,]),
('3ICQ_TP0_backrub', [1,]),
('3M1I_TP2_backrub', [1,]),
('3W3Z_TP0_backrub', [1,])
]
}
# # everything else
# 5 : [
# ('1A2K_TP0_backrub', [1,]),
# ('1WA5_TP1_backrub', [1,]),
# ('3A6P_TP0_backrub', [1,]),
# ('3EA5_TP0_backrub', [1,]),
# ('3M1I_TP1_backrub', [1,]),
# ('3WYF_TP1_backrub', [1,]),
# ('4OL0_TP0_backrub', [1,])
# ]
# }

# # macrostates - dict used for organizing other data structures
# macrostates = {
# # GEF-complexed, no nucleotide - one structure: 1I2M
# 1 : [
# ('1I2M_TP0_backrub', [1,2,3,4,5,6])
# ]
# }

target_pos = [101,102,103,105,106,107,108,109,122,123,124,125,126,127]
seq_to_arr_idx = {}
for idx,pos in enumerate(target_pos):
    seq_to_arr_idx[pos] = idx

#print seq_to_arr_idx

# turn the index file into a dictionary for much faster lookups
seq_dict = {}
for index, row in seq_index.iterrows():
    seq_dict[(row['pdb'], row['pdb_res_num'])] = row['yeast_seq_num']

# given a structure_id and a res number in pdb numbering, return the yeast seq number and wt aa id
def lookup_res(structure_id, pdb_num):
    try:
        yeast_num = seq_dict[(structure_id.lower(), pdb_num)]
        wt_aa_id = yeast_gsp1_seq[yeast_num-1]
        return yeast_num, wt_aa_id
    except KeyError:
        print('structure id not found in index file')
        raise

# parse and organize all the score files
score_data = {}
for key in macrostates:
    for substate in macrostates[key]:
        for i in substate[1]:
            substate_name = '%s_%s' % (substate[0], i)
            if os.path.isfile(os.path.join(scorefile_dir, '%s.sc' % substate_name)):
                # store a tuple for each substate, where the first value is the macrostate id,
                # the second is a list of the backrub ensemble ids that went into the substate,
                # and third is the dataframe containing scores for the substate
                with open(os.path.join(scorefile_dir, '%s.sc' % substate_name), 'r') as data_file:
                    data_lines = data_file.readlines()
                data_lines = [i.strip('\n').split() for i in data_lines[2:]]
                ensemble_ids = set()
                for line in data_lines:
                    iden = line[-1].split('_')
                    ensemble_ids.add(int(iden[2][7:]))
                score_data[substate_name] = (key, sorted(ensemble_ids), data_lines)

# for key in sorted(score_data):
#     print key, score_data[key][0], score_data[key][1], len(score_data[key][2])
# print len(score_data)

# organize all the scores for a specified macrostate into an array for processing
def build_mstate_array(macrostate_id):
    # build dict for indexing data array
    datalabels = {}
    datalabel_names = []
    for substate in score_data:
        if score_data[substate][0] == macrostate_id:
            datalabel_names.append(substate)
    for idx,substate in enumerate(sorted(datalabel_names)):
        datalabels[substate] = (idx, {})
        for jdx,ustate_id in enumerate(score_data[substate][1]):
            datalabels[substate][1][ustate_id] = jdx

    # for datalabel in sorted(datalabels):
    #     print datalabel, datalabels[datalabel][0]
    #     for us_id in sorted(datalabels[datalabel][1]):
    #         print us_id, datalabels[datalabel][1][us_id]

    # build an empty array to hold data
    # 1st index is substate id
    # 2nd index is microstate id
    # 3nd index is mutation
    # 4th index is residue number
    datarray = np.zeros((len(datalabels), 10, 20, 14))
    # fill the array up, iterating over each substate
    for ss_id in datalabels:
        # get corresponding data
        print('processing %s scores' % ss_id)
        for line in score_data[ss_id][2]:
            iden = line[-1].split('_')
            rel_seq_idx, nataa = lookup_res('%s_%s' % (iden[0], iden[1]), int(iden[4][:-1]))
            #print rel_seq_idx, nataa
            print('%s_%s' % (iden[0], iden[1]), iden[4][:-1], iden[4][-1], line[1], rel_seq_idx, nataa)
            # print iden
            # insert score into the array
            # print datalabels[ss_id][0]
            # print datalabels[ss_id][1][int(iden[2][7:])]
            # print aminotonumber[iden[4][-1]]-1
            # print seq_to_arr_idx[int(rel_seq_idx)]
            datarray[datalabels[ss_id][0],
                     datalabels[ss_id][1][int(iden[2][7:])],
                     aminotonumber[iden[4][-1]]-1,
                     seq_to_arr_idx[int(rel_seq_idx)]] = float(line[1])
    return (datarray, datalabels)

for state in macrostates:
    #print build_mstate_array(state)
    state_data_array = build_mstate_array(state)
    pickle.dump(state_data_array, open('state_%s_data_v3.pkl' % str(state), 'wb'))
