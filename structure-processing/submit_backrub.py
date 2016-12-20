#!/usr/bin/python
#$ -S /usr/bin/python
#$ -cwd
#$ -r yes
#$ -j y
#$ -l h_rt=48:00:00
#$ -t 1-200
#$ -l arch=linux-x64
#$ -l mem_free=2G

import socket
import sys
import os
import subprocess
import shlex

rosetta_backrub = "/Users/anatale/Rosetta/main/source/bin/backrub.default.macosclangrelease"
rosetta_db_dir = "/Users/anatale/Rosetta/main/database"

input_pdbs_dir = "/Users/anatale/Documents/school/UCSF/Kortemme_lab/rosetta_staging/testing"
listfile = os.path.join(input_pdbs_dir, 'params_list.txt')

ntrials=100

sge_task_id = 1
if os.environ.has_key("SGE_TASK_ID"):
	sge_task_id = os.environ["SGE_TASK_ID"]

output_dir = "/Users/anatale/Documents/school/UCSF/Kortemme_lab/rosetta_staging/testing/output"
outfile_name = os.path.join(output_dir, 'backrub_log_%s.log' % str(sge_task_id))

def generate_backrub_cmd(pdb_file, params_file, pivot_residues):
    if params_file == 'NA':
        backrub_cmd = shlex.split(rosetta_backrub+
        ' -database '+rosetta_db_dir+
        ' -s '+os.path.join(input_pdbs_dir, pdb_file)+
        ' -ignore_waters -ex1 -ex2 -extrachi_cutoff 0 '+
        ' -backrub:ntrials '+str(ntrials)+
        ' -nstruct 1 -mc_kt 0.9 -backrub:initial_pack '+
        ' -backrub:pivot_residues '+pivot_residues+
        ' -out:no_nstruct_label '+
        ' -out:suffix '+str(sge_task_id))
    else:
        backrub_cmd = shlex.split(rosetta_backrub+
        ' -database '+rosetta_db_dir+
        ' -s '+os.path.join(input_pdbs_dir, pdb_file)+
        ' -extra_res_fa '+os.path.join(input_pdbs_dir, params_file)+
        ' -ignore_waters -ex1 -ex2 -extrachi_cutoff 0 '+
        ' -backrub:ntrials '+str(ntrials)+
        ' -nstruct 1 -mc_kt 0.9 -backrub:initial_pack '+
        ' -backrub:pivot_residues '+pivot_residues+
        ' -out:no_nstruct_label '+
        ' -out:suffix '+str(sge_task_id))
    return backrub_cmd

with open(outfile_name, 'w') as outfile:
    outfile.write("Python: %s\n" % sys.version)
    outfile.write("Host: %s\n" % socket.gethostname())

with open(listfile, 'r') as params_list:
    pdb_params = params_list.readlines()
pdb_params = [i.strip('\n') for i in pdb_params]

for line in pdb_params[1:]:
    params = line.split()
    backrub_outdir = os.path.join(output_dir, params[0][:-4])
    try:
        os.makedirs(backrub_outdir)
    except OSError:
        if not os.path.isdir(backrub_outdir):
            raise
    os.chdir(backrub_outdir)
    backrub_cmd = generate_backrub_cmd(params[0], params[1], \
      ' '.join([str(res) for res in params[2:]]))
    with open(outfile_name, 'a+') as outfile:
        outfile.write('\n*************************************************\n')
        outfile.write('Starting backrub job for pdb file: %s' % params[0])
        outfile.write('\n*************************************************\n\n')
    with open(outfile_name, 'a+') as outfile:
        process = subprocess.Popen(backrub_cmd, \
                                   stdout=outfile, \
                                   stderr=subprocess.STDOUT, \
                                   close_fds = True)
        returncode = process.wait()
    # cleanup
    for filename in os.listdir(os.getcwd()):
        if not filename.endswith('_last.pdb'):
            os.remove(filename)
    os.chdir(output_dir)
