#!/usr/bin/env python2

import random
import os
import sys
import copy

ensemble_size = 200
x = range(1, ensemble_size+1)

# _dir -> exact location to find or write files
# _prefix -> base dir path, needs a specified subdir
# output:
output_dir = '/Users/anatale/school/UCSF/Kortemme_lab/rosetta_staging/state_1a2k/state_lists3'
final_list_path = '/netapp/home/anatale/fixbb_run3/state_lists3'
blacklists_dir = '/Users/anatale/school/UCSF/Kortemme_lab/rosetta_staging/blacklist'
# input:
pdb_path_prefix = '/netapp/home/anatale/backrub_out'
paramsfiles_dir = '/netapp/home/anatale/fixbb_run3/params'
resfiles_prefix = '/netapp/home/anatale/fixbb_run3/resfiles3'
# make output directory
try:
    os.makedirs(output_dir)
except OSError:
    if not os.path.isdir(output_dir):
        raise

# macrostates
# macrostates = {
# # GEF-complexed, no nucleotide - one structure: 1I2M
# 1 : [
# ('1I2M_TP0_backrub', 'NA', 40, 20, '1I2M_TP0_backrub_blacklist'), # already have 60
# ],
# # GAP-complexed, GTP-analogue nucleotide - one structure: 1K5D
# 2 : [
# ('1K5D_TP0_backrub', '1K5D_TP0_backrub.params', 40, 20, '1K5D_TP0_backrub_blacklist'), # already have 60
# ],
# # GAP and BP1-complexed, GTP-analogue nucleotide - one structure: 1K5D
# 3 : [
# ('1K5D_ALL_backrub', '1K5D_ALL_backrub.params', 100, 20, None),
# ],
# # uncomplexed, GDP nucleotide - one structure: 3GJ0
# 4 : [
# ('3GJ0_TP0_backrub', '3GJ0_TP0_backrub.params', 80, 20, '3GJ0_TP0_backrub_blacklist'), # already have 20
# ],
# # karyopherin complexed, GTP nucleotide - four structures: 1QBK, 2BKU, 3W3Z, 4OL0
# 5 : [
# ('1QBK_TP0_backrub', '1QBK_TP0_backrub.params', 15, 15, '1QBK_TP0_backrub_blacklist'), # already have 10
# ('2BKU_TP0_backrub', '2BKU_TP0_backrub.params', 15, 15, '2BKU_TP0_backrub_blacklist'), # already have 10
# ('3W3Z_TP0_backrub', '3W3Z_TP0_backrub.params', 15, 15, '3W3Z_TP0_backrub_blacklist'), # already have 10
# ('4OL0_TP0_backrub', '4OL0_TP0_backrub.params', 25, 25, None)
# ],
# # modified structures with no bp - combine all the above states
# 6 : [
# ('1I2M_NON_backrub', 'NA', 25, 25, None),
# ('1K5D_NON_backrub', '1K5D_NON_backrub.params', 25, 25, None),
# ('1QBK_NON_backrub', '1QBK_NON_backrub.params', 25, 25, None),
# ('3GJ0_TP0_backrub0', '3GJ0_TP0_backrub.params', 25, 25, None) # let this repick any state from the ensemble
# ]
# }

# # macrostates for testing
# macrostates = {
# 1 : [
# ('2BKU_TP0_backrub', '2BKU_TP0_backrub.params', 60),
# ]
# }

# macrostates
macrostates = {
1 : [
('1A2K_TP0_backrub', '1A2K_TP0_backrub.params', 100, 10, None), # already have 60
],
2 : [
('1A2K_NON_backrub', '1A2K_NON_backrub.params', 25, 25, None), # already have 60
]
}
# for each of the structures defined above, write out lists of randomly
# selected backrub pdbs which will be processed as
def write_lists(options, list_params):
    # make a list of random intergers to pull from (set up above)
    y = copy.copy(options)
    structure_name, paramsfile, count, entries_per_list, blacklist = list_params
    if blacklist:
        with open(os.path.join(blacklists_dir, blacklist), 'r') as blacklist_file:
            reader = blacklist_file.readline()
        for n in reader.split():
            if int(n) in y:
                y.remove(int(n))
    random.shuffle(y)
    list_number = 1
    list_of_listfiles = []
    while entries_per_list * (list_number - 1) < count:
        # name the list with the structure it comes from and a number tag
        list_file_name = '%s_%s.lst' % (structure_name, str(list_number))
        with open(os.path.join(output_dir, list_file_name), 'w') as pdb_list:
            # take chunks out of the list of random ints to select pdb files
            for i in y[(list_number-1)*entries_per_list:list_number*entries_per_list]:
                target_pdb_path = os.path.join(pdb_path_prefix,
                                               '%s' % structure_name,
                                               '%s%s_last.pdb' % (structure_name, i))
                # # make sure target pdb exists, then write it to the list
                # if os.path.isfile(target_pdb_path):
                #     pdb_list.write('%s\n' % target_pdb_path)
                # else:
                #     print 'target pdb does not exist!'
                #     print target_pdb_path
                #     sys.exit()
                pdb_list.write('%s\n' % target_pdb_path)
        # save the names of all the lists made for this structure
        list_of_listfiles.append(os.path.join(final_list_path, list_file_name))
        list_number += 1
    return list_of_listfiles

# read the blacklist, then return the randomized number list without those values
def test_blacklist(options, blacklist):
    print '-------------------------------------------------'
    z = copy.copy(options)
    print len(z), '\n'
    with open(os.path.join(blacklists_dir, blacklist), 'r') as blacklist_file:
        reader = blacklist_file.readline()
    print len(reader.split()), '\n'
    for n in reader.split():
        #print n
        if int(n) in z:
            z.remove(int(n))
    print len(z)

def process_states():
    job_id = 1
    # make a list of jobs to be read by the final script
    with open(os.path.join(output_dir, 'fixbb_jobs.txt'), 'w') as jobsfile:
        jobsfile.write('job_id pdb_list paramsfile resfile_list\n')
        # write lists of pdbs
        for state in macrostates:
            for list_params in macrostates[state]:
                list_of_listfiles = write_lists(x, list_params)
                # write all the job parameters into the jobs list file
                if list_params[1] == 'NA':
                    paramsfile_path = 'NA'
                else:
                    paramsfile_path = os.path.join(paramsfiles_dir, list_params[1])
                resfile_list_path = os.path.join(resfiles_prefix, list_params[0], '%s_resfiles.lst' % list_params[0])
                for list_name in list_of_listfiles:
                    jobsfile.write('%s %s %s %s\n' % (str(job_id), list_name, paramsfile_path, resfile_list_path))
                    job_id += 1

process_states()
# for state in macrostates:
#     for substate in macrostates[state]:
#         blacklist = substate[3]
#         if blacklist:
#             test_blacklist(x, blacklist)
