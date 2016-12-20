#/usr/bin/env python2
import sys
import pandas as pd
import os
import numpy as np
import pickle

seq_index = pd.read_table('/Users/anatale/Documents/school/UCSF/Kortemme_lab/code/structure-processing/aux_files/Gsp1_index_all.txt')

yeast_gsp1_seq = '\
MSAPAANGEVPTFKLVLVGDGGTGKTTFVKRHLTGEFEKKYIATIGVEVHPLSFYTNFGE\
IKFDVWDTAGQEKFGGLRDGYYINAQCAIIMFDVTSRITYKNVPNWHRDLVRVCENIPIV\
LCGNKVDVKERKVKAKTITFHRKKNLQYYDISAKSNYNFEKPFLWLARKLAGNPQLEFVA\
SPALAPPEVQVDEQLMQQYQQEMEQATALPLPDEDDADL'

# given a structure_id and a res number in pdb numbering, return the yeast seq number and wt aa id
def lookup_res(structure_id, pdb_num):
    if structure_id.lower() in seq_index['pdb'].tolist():
        for index, row in seq_index.iterrows():
            if row['pdb'] == structure_id.lower() and row['pdb_res_num'] == pdb_num:
                yeast_num = int(row['yeast_seq_num'])
    else:
        print 'structure id not found in index file'
        sys.exit()
    wt_aa_id = yeast_gsp1_seq[yeast_num-1]
    return yeast_num, wt_aa_id

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
],
# everything else
5 : [
('1A2K_TP0_backrub', [1,]),
('1WA5_TP1_backrub', [1,]),
('3A6P_TP0_backrub', [1,]),
('3EA5_TP0_backrub', [1,]),
('3M1I_TP1_backrub', [1,]),
('3WYF_TP1_backrub', [1,]),
('4OL0_TP0_backrub', [1,])
]
}

# load all score files into dataframes with some id info
score_dataframes = {}
for key in macrostates:
    for substate in macrostates[key]:
        for i in substate[1]:
            substate_name = '%s_%s' % (substate[0], i)
            if os.path.isfile(os.path.join(os.getcwd(), '%s.sc' % substate_name)):
                # store a tuple for each substate, where the first value is the macrostate id,
                # the second is a list of the backrub ensemble ids that went into the substate,
                # and third is the dataframe containing scores for the substate
                score_dataframes[substate_name] = (key, [], pd.read_csv('%s.sc' % substate_name, delim_whitespace=True, header=1))
            ensemble_ids = set()
            for idx, row in score_dataframes[substate_name][2].iterrows():
                iden = row['description'].split('_')
                ensemble_ids.add(int(iden[2][7:]))
            score_dataframes[substate_name][1].extend(sorted(ensemble_ids))

for key in sorted(score_dataframes):
    print key, score_dataframes[key][0], score_dataframes[key][1], score_dataframes[key][2].size
print len(score_dataframes)

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

# organize all the scores for a specified macrostate into an array for processing
def build_mstate_array(macrostate_id):
    # build dict for indexing data array
    datalabels = {}
    datalabel_names = []
    for substate in score_dataframes:
        if score_dataframes[substate][0] == macrostate_id:
            datalabel_names.append(substate)
    for idx,substate in enumerate(sorted(datalabel_names)):
        datalabels[substate] = (idx, {})
        for jdx,ustate_id in enumerate(score_dataframes[substate][1]):
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
    datarray = np.zeros((len(datalabels), 10, 20, 6))
    # fill the array up, iterating over each substate
    for ms_id in datalabels:
        # get corresponding dataframe
        print 'processing %s scores' % ms_id
        for i,row in score_dataframes[ms_id][2].iterrows():
            iden = row['description'].split('_')
            #print '%s_%s' % (iden[0], iden[1]), iden[4][:-1]
            rel_seq_idx, nataa = lookup_res('%s_%s' % (iden[0], iden[1]), int(iden[4][:-1]))
            # insert score into the array
            datarray[datalabels[ms_id][0],
                     datalabels[ms_id][1][int(iden[2][7:])],
                     aminotonumber[iden[4][3]]-1,
                     rel_seq_idx-122] = float(row['total_score'])
    return (datarray, datalabels)

for state in macrostates:
    state_data_array = build_mstate_array(state)
    pickle.dump(state_data_array, open('state_%s_data.pkl' % str(state), 'wb'))
