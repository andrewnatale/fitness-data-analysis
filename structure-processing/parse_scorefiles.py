import sys
import pandas as pd
import os
import numpy as np
import pickle
import csv

scorefile_dir = '/Users/anatale/school/UCSF/Kortemme_lab/data/combo_scores'

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
('1K5D_TP0_backrub', [1,2,3,4,5,6,7,8])
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
('3GJ0_TP0_backrub', [7,])
]
}

target_pos = [101,102,103,105,106,107,108,109,122,123,124,125,126,127]
seq_to_arr_idx = {}
for idx,pos in enumerate(target_pos):
    seq_to_arr_idx[pos] = idx

# turn the index file into a dictionary for much faster lookups
seq_dict = {}
for index, row in seq_index.iterrows():
    seq_dict[(row['pdb'][:4].upper(), row['pdb_res_num'])] = row['yeast_seq_num']

# given a structure_id and a res number in pdb numbering, return the yeast seq number and wt aa id
def lookup_res(structure_id, pdb_num):
    try:
        yeast_num = seq_dict[(structure_id, pdb_num)]
        #wt_aa_id = yeast_gsp1_seq[yeast_num-1]
        return yeast_num
    except KeyError:
        print('structure id not found in index file')
        raise


class MacroStateData(object):

    def __init__(self, state_identifiers):
        # generate scorefile names to load
        self.substates = []
        for ss_name in state_identifiers:
            for ss_id in ss_name[1]:
                self.substates.append('%s_%d' % (ss_name[0], ss_id))
        #print(self.substates)

    def load_scorefiles(self):
        self.ss_data_list = []
        for ss_name in self.substates:
            scorefile_name = os.path.join(scorefile_dir, '%s.sc' % ss_name)
            if os.path.isfile(scorefile_name):
                self.ss_data_list.append(SubStateData(ss_name, scorefile_name))
            else:
                print('score file not found, exiting!')
                exit(0)
        #print(self.ss_data_list)

    def generate_array(self):
        us_count = 0
        self.us_idx_dict = {}
        for ss in self.ss_data_list:
            for score in ss.scorelist:
                if not (ss.name, score[1]) in self.us_idx_dict:
                    self.us_idx_dict[(ss.name, score[1])] = us_count
                    us_count += 1
        self.state_array = np.zeros((us_count, 20, len(target_pos)))
        self.array_labels = []
        for ss in self.ss_data_list:
            for score in ss.scorelist:
                us_idx = self.us_idx_dict[(ss.name, score[1])]
                mut_idx = aminotonumber[score[3]] - 1
                pos_idx = seq_to_arr_idx[lookup_res(ss.name[:4], int(score[2]))]
                self.state_array[us_idx, mut_idx, pos_idx] = float(score[0])
        print(self.state_array.size)

    def write_pkl(self, outname):
        pickle.dump((self.state_array, self.us_idx_dict), open(outname, 'wb'))

class SubStateData(object):

    def __init__(self, ss_name, scorefile_name):
        self.name = ss_name
        self.backrub_ids = set()
        self.scorelist = []
        with open(scorefile_name, 'r') as scorefile:
            reader = scorefile.readlines()
        for line in [i.strip('\n').split() for i in reader[2:]]:
            total_score = line[1]
            score_id = line[-1].split('_')
            backrub_id = score_id[2][7:]
            position_id = score_id[-1][:-1]
            mutation_id = score_id[-1][-1]
            #print(backrub_id, position_id, mutation_id, total_score)
            self.backrub_ids.add(backrub_id)
            self.scorelist.append((total_score, backrub_id, position_id, mutation_id))
        # print(self.name)
        # print(sorted(self.backrub_ids))
        # print(len(self.scorelist))

all_data = []
for state in macrostates:
    all_data.append(MacroStateData(macrostates[state]))
for idx,obj in enumerate(all_data):
    obj.load_scorefiles()
    obj.generate_array()
    obj.write_pkl('state_%s_scores.pkl' % str(idx+1))
