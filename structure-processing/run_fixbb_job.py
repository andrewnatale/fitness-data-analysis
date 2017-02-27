#!/usr/bin/env
#$ -S /usr/bin/python
#$ -cwd
#$ -r yes
#$ -j y
#$ -l h_rt=12:00:00
#$ -t 1-280
#$ -l arch=linux-x64
#$ -l mem_free=2G

# This script takes a list of pdb files, a params file, and a numbered list
# of res files as input. It will run rosetta on the whole list of pdbs using
# the res file in the list corresponding to the sge_task_id variable

import socket
import sys
import os
import subprocess
import shlex

script, job_pdb_list, job_paramsfile, job_resfile_list = sys.argv

rosetta_fixbb = "/netapp/home/gkreder/Rosetta_new/source/bin/fixbb.default.linuxgccrelease"
rosetta_db_dir = "/netapp/home/gkreder/Rosetta_new/database"

resfiles_prefix = '/netapp/home/anatale/fixbb_run3/resfiles'

sge_task_id = 1
if os.environ.has_key("SGE_TASK_ID"):
	sge_task_id = os.environ["SGE_TASK_ID"]

job_struct = job_pdb_list.split('/')[-1][:-4]

# find the desired resfile according to the sge_task_id variable
with open(job_resfile_list, 'r') as infile2:
    resfiles = infile2.readlines()
resfiles = [i.strip('\n') for i in resfiles]

task_resfile = None
mutation_tag = None
for resfile in resfiles:
    tmp_resfile_id, tmp_resfile = resfile.split()
    if str(tmp_resfile_id) == str(sge_task_id):
        task_resfile = os.path.join(resfiles_prefix, job_struct[:-2], tmp_resfile)
        mutation_tag = tmp_resfile.split('_')[3].split('.')[0]
        break

output_prefix = '/netapp/home/anatale/fixbb_run3/output'
output_suffix = job_struct

# make output directory
try:
    os.makedirs(os.path.join(output_prefix, output_suffix, mutation_tag))
except OSError:
    if not os.path.isdir(os.path.join(output_prefix, output_suffix, mutation_tag)):
        raise

def generate_fixbb_cmd(pdb_listfile, resfile, params_path):
    if params_path == 'NA':
        fixbb_cmd = shlex.split(rosetta_fixbb+
        ' -database '+rosetta_db_dir+
        ' -l '+pdb_listfile+
        ' -resfile '+resfile+
        ' -ignore_waters -ex1 -ex2 -extrachi_cutoff 0'+
        ' -nstruct 1'+
        ' -linmem_ig 10'+
        ' -ignore_unrecognized_res'+
        ' -minimize_sidechains'+
        ' -multi_cool_annealer 10'+
        ' -out:no_nstruct_label')
    else:
        fixbb_cmd = shlex.split(rosetta_fixbb+
        ' -database '+rosetta_db_dir+
        ' -l '+pdb_listfile+
        ' -resfile '+resfile+
        ' -extra_res_fa '+params_path+
        ' -ignore_waters -ex1 -ex2 -extrachi_cutoff 0'+
        ' -nstruct 1'+
        ' -linmem_ig 10'+
        ' -ignore_unrecognized_res'+
        ' -minimize_sidechains'+
        ' -multi_cool_annealer 10'+
        ' -out:no_nstruct_label')
    return fixbb_cmd

logfile_name = os.path.join(output_prefix, output_suffix, 'fixbb_%s-%s.log' % (output_suffix, mutation_tag))

# chdir into output file directory
os.chdir(os.path.join(output_prefix, output_suffix, mutation_tag))

with open(logfile_name, 'w') as outfile:
    outfile.write("Python: %s\n" % sys.version)
    outfile.write("Host: %s\n" % socket.gethostname())

fixbb_cmd = generate_fixbb_cmd(job_pdb_list, task_resfile, job_paramsfile)
#print fixbb_cmd

with open(logfile_name, 'a+') as outfile:
    outfile.write('\n*************************************************\n')
    outfile.write('Starting fixbb job for pdb file list: %s\n' % job_pdb_list)
    outfile.write('Using resfile: %s' % task_resfile)
    outfile.write('\n*************************************************\n\n')
with open(logfile_name, 'a+') as outfile:
    process = subprocess.Popen(fixbb_cmd, \
                               stdout=outfile, \
                               stderr=subprocess.STDOUT, \
                               close_fds = True)
    returncode = process.wait()

for pdb_file in os.listdir(os.getcwd()):
    if pdb_file.endswith('.pdb'):
        os.rename(pdb_file, '%s_fixbb_%s.pdb' % (pdb_file[:-9], mutation_tag))
