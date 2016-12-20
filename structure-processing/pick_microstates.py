#!/usr/bin/env python2

import random
import os
import sys

ensemble_size = 200
x = range(1, ensemble_size+1)

# _dir -> exact location to find or write files
# _prefix -> base dir path, needs a specified subdir
# output:
output_dir = '/netapp/home/anatale/fixbb/state_lists'
# input:
pdb_path_prefix = '/netapp/home/anatale/backrub_out'
paramsfiles_dir = '/netapp/home/anatale/fixbb/params'
resfiles_prefix = '/netapp/home/anatale/fixbb/resfiles'
# make output directory
try:
    os.makedirs(output_dir)
except OSError:
    if not os.path.isdir(output_dir):
        raise

# macrostates
macrostates = {
# GEF-complexed, no nucleotide - one structure: 1I2M
1 : [
('1I2M_TP0_backrub', 'NA', 60)
],
# GAP-complexed, GTP-analogue nucleotide - one structure: 1K5D
2 : [
('1K5D_TP0_backrub', '1K5D_TP0_backrub.params', 60)
],
# uncomplexed, GDP nucleotide - three structures: 3GJ0, 5CIQ, 5CIT
3 : [
('3GJ0_TP0_backrub', '3GJ0_TP0_backrub.params', 20),
('5CIQ_TP0_backrub', '5CIQ_TP0_backrub.params', 20),
('5CIT_TP0_backrub', '5CIT_TP0_backrub.params', 20)
],
# karyopherin complexed, GTP nucleotide - six structures: 1QBK, 1WA5_TP2, 2BKU, 3ICQ, 3M1I_TP2, 3W3Z
4 : [
('1QBK_TP0_backrub', '1QBK_TP0_backrub.params', 10),
('1WA5_TP2_backrub', '1WA5_TP2_backrub.params', 10),
('2BKU_TP0_backrub', '2BKU_TP0_backrub.params', 10),
('3ICQ_TP0_backrub', '3ICQ_TP0_backrub.params', 10),
('3M1I_TP2_backrub', '3M1I_TP2_backrub.params', 10),
('3W3Z_TP0_backrub', '3W3Z_TP0_backrub.params', 10)
]
}
# # everything else
# 5 : [
# ('1A2K_TP0_backrub', '1A2K_TP0_backrub.params', 10),
# ('1WA5_TP1_backrub', '1WA5_TP1_backrub.params', 10),
# ('3A6P_TP0_backrub', '3A6P_TP0_backrub.params', 10),
# ('3EA5_TP0_backrub', '3EA5_TP0_backrub.params', 10),
# ('3M1I_TP1_backrub', '3M1I_TP1_backrub.params', 10),
# ('3WYF_TP1_backrub', '3WYF_TP1_backrub.params', 10),
# ('4OL0_TP0_backrub', '4OL0_TP0_backrub.params', 10)
# ]
# }

# # macrostates for testing
# macrostates = {
# 1 : [
# ('2BKU_TP0_backrub', '2BKU_TP0_backrub.params', 60),
# ]
# }

# for each of the structures defined above, write out lists of randomly
# selected backrub pdbs which will be processed as
def write_lists(structure_name, count, options):
    # make a list of random intergers to pull from (set up above)
    y = options
    random.shuffle(y)
    # how many structures per list?
    entries_per_list = 10
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
                # make sure target pdb exists, then write it to the list
                if os.path.isfile(target_pdb_path):
                    pdb_list.write('%s\n' % target_pdb_path)
                else:
                    print 'target pdb does not exist!'
                    print target_pdb_path
                    sys.exit()

        # save the names of all the lists made for this structure
        list_of_listfiles.append(os.path.join(output_dir, list_file_name))
        list_number += 1
    return list_of_listfiles

job_id = 1
# make a list of jobs to be read by the final script
with open(os.path.join(output_dir, 'fixbb_jobs.txt'), 'w') as jobsfile:
    jobsfile.write('job_id pdb_list paramsfile resfile_list\n')
    # write lists of pdbs
    for state in macrostates:
        for structure in macrostates[state]:
            list_of_listfiles = write_lists(structure[0], structure[2], x)
            # write all the job parameters into the jobs list file
            if structure[1] == 'NA':
                paramsfile_path = 'NA'
            else:
                paramsfile_path = os.path.join(paramsfiles_dir, structure[1])
            resfile_list_path = os.path.join(resfiles_prefix, structure[0], '%s_resfiles.lst' % structure[0])
            for list_name in list_of_listfiles:
                jobsfile.write('%s %s %s %s\n' % (str(job_id), list_name, paramsfile_path, resfile_list_path))
                job_id += 1
