#/usr/bin/env python2
import sys
import pandas as pd
import os
import pickle

seq_index = pd.read_table('/Users/anatale/Documents/school/UCSF/Kortemme_lab/code/structure-processing/Gsp1_index_all.txt')

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

backrub_scores = {}
for filename in os.listdir(os.getcwd()):
    if filename.endswith('.sc'):
        df = pd.read_csv(filename, delim_whitespace=True, header=1)
        for idx, row in df.iterrows():
            iden = row['description'].split('_')
            #print '%s_%s' % (iden[0], iden[1]), iden[2][7:], row['total_score']
            backrub_scores[('%s_%s' % (iden[0], iden[1]), iden[2][7:])] = row['total_score']

# for iden in backrub_scores:
#     print iden, backrub_scores[iden]

pickle.dump(backrub_scores, open('backrub_scores.pkl', 'wb'))
