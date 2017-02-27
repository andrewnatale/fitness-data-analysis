#!/usr/bin/python
#$ -S /usr/bin/python
#$ -cwd
#$ -r yes
#$ -j y
#$ -l h_rt=48:00:00
#$ -t 1-21
#$ -l arch=linux-x64
#$ -l mem_free=2G

# for each of the fixbb jobs defined in the job list, write a list of all
# the output pdbs for scoring (1200 pdbs in each list), then score them all together

import socket
import sys
import os
import subprocess
import shlex

rosetta_score_jd2 = "/netapp/home/gkreder/Rosetta_new/source/bin/score_jd2.default.macosclangrelease"
rosetta_db_dir = "/netapp/home/gkreder/Rosetta_new/database"

job_list_path = '/netapp/home/anatale/fixbb_run3/state_lists/fixbb_jobs.txt'
pdb_search_path = '/netapp/home/anatale/fixbb_run3/output'
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

with open(job_list_path, 'r') as infile:
    jobs = infile.readlines()
jobs = [i.strip('\n') for i in jobs]

# variables to unpack from job list
job_id = None
job_pdb_list = None
job_paramsfile = None
job_resfile_list = None

job_name = None

# unpack all the paths and variables needed for this job
for job in jobs[1:]:
    if str(sge_task_id) == str(job.split()[0]):
        job_id, job_pdb_list, job_paramsfile, job_resfile_list = job.split()
        if job_id and job_pdb_list and job_paramsfile and job_resfile_list:
            job_name = job_pdb_list.split('/')[-1][:-4]
        else:
            print 'cannot start job, id not found!'
            sys.exit()

# first write a list of all the pdbs in a format that can be passed to Rosetta
score_list_path = os.path.join(output_prefix, '%s_fixbb_pdbs.lst' % job_name)
with open(score_list_path, 'w') as score_list:
    for dirpath,dirnames,filenames in os.walk(os.path.join(pdb_search_path, job_name)):
        for filename in filenames:
            if filename.endswith('.pdb'):
                score_list.write('%s\n' % os.path.join(dirpath, filename))

logfile_name = os.path.join(output_prefix, '%s-scores.log' % job_name)

# use list to run Rosetta
score_cmd = generate_score_cmd(score_list_path, job_paramsfile, job_name)
print score_cmd

with open(logfile_name, 'w') as outfile:
    outfile.write("Python: %s\n" % sys.version)
    outfile.write("Host: %s\n" % socket.gethostname())

with open(logfile_name, 'a+') as outfile:
    outfile.write('\n*************************************************\n')
    outfile.write('Starting scoring job for file list: %s\n' % job_pdb_list)
    outfile.write('\n*************************************************\n\n')
with open(logfile_name, 'a+') as outfile:
    process = subprocess.Popen(score_cmd, \
                               stdout=outfile, \
                               stderr=subprocess.STDOUT, \
                               close_fds = True)
    returncode = process.wait()
