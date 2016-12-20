#!/usr/bin/python
#$ -S /usr/bin/python
#$ -cwd
#$ -r yes
#$ -j y
#$ -l h_rt=48:00:00
#$ -t 1-18
#$ -l arch=linux-x64
#$ -l mem_free=2G

# score each backrub structure

import socket
import sys
import os
import subprocess
import shlex

rosetta_score_jd2 = "/Users/anatale/Rosetta/main/source/bin/score_jd2.default.macosclangrelease"
rosetta_db_dir = "/Users/anatale/Rosetta/main/database"

#job_list_path = '/Users/anatale/Documents/school/UCSF/Kortemme_lab/rosetta_staging/state_lists/fixbb_jobs.txt'
pdb_search_path = '/Users/anatale/Documents/school/UCSF/Kortemme_lab/data'
params_search_path = '/Users/anatale/Documents/school/UCSF/Kortemme_lab/structures/backrub_ready'
output_prefix = os.getcwd()

sge_task_id = 1
if os.environ.has_key("SGE_TASK_ID"):
	sge_task_id = os.environ["SGE_TASK_ID"]

def generate_score_cmd(pdb_listfile, params_path, output_name):
    if params_path == 'NA':
        score_cmd = shlex.split(rosetta_score_jd2+\
        ' -database '+rosetta_db_dir+\
        ' -in:file:l '+pdb_listfile+\
        ' -ignore_unrecognized_res '+\
        ' -out:no_nstruct_label'+\
        ' -out:file:score_only'+\
        ' -out:file:scorefile %s.sc' % output_name)
    else:
        score_cmd = shlex.split(rosetta_score_jd2+\
        ' -database '+rosetta_db_dir+\
        ' -in:file:l '+pdb_listfile+\
        ' -extra_res_fa '+params_path+\
        ' -ignore_unrecognized_res'+\
        ' -out:no_nstruct_label'+\
        ' -out:file:score_only'+\
        ' -out:file:scorefile %s.sc' % output_name)
    return score_cmd

structure_list = sorted([dirname for dirname in os.listdir(pdb_search_path) if dirname.endswith('_backrub')])

job_name = structure_list[int(sge_task_id)-1]

# find params file for this structure
job_paramsfile = 'NA'
for filename in [filename for filename in os.listdir(params_search_path) if filename.endswith('.params')]:
    if filename.split('.')[0] == job_name:
        job_paramsfile = os.path.join(params_search_path, filename)

# first write a list of all the pdbs in a format that can be passed to Rosetta
score_list_path = os.path.join(output_prefix, '%s_all_pdbs.lst' % job_name)
with open(score_list_path, 'w') as score_list:
    for dirpath,dirnames,filenames in os.walk(os.path.join(pdb_search_path, job_name)):
        for filename in filenames:
            if filename.endswith('.pdb'):
                score_list.write('%s\n' % os.path.join(dirpath, filename))

logfile_name = os.path.join(output_prefix, '%s-scores.log' % job_name)

# use list to run Rosetta
score_cmd = generate_score_cmd(score_list_path, job_paramsfile, job_name)
#print score_cmd

with open(logfile_name, 'w') as outfile:
    outfile.write("Python: %s\n" % sys.version)
    outfile.write("Host: %s\n" % socket.gethostname())

with open(logfile_name, 'a+') as outfile:
    outfile.write('\n*************************************************\n')
    outfile.write('Starting scoring job for file list: %s\n' % job_name)
    outfile.write('\n*************************************************\n\n')
with open(logfile_name, 'a+') as outfile:
    process = subprocess.Popen(score_cmd, \
                               stdout=outfile, \
                               stderr=subprocess.STDOUT, \
                               close_fds = True)
    returncode = process.wait()
