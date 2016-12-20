#!/usr/bin/env python2

import socket
import sys
import os
import subprocess
import shlex

fixbb_script = sys.argv[1]
job_list = sys.argv[2]

# find the arguments needed to run jobs
with open(job_list, 'r') as infile:
    jobs = infile.readlines()
jobs = [i.strip('\n') for i in jobs]

job_id = None
job_pdb_list = None
job_paramsfile = None
job_resfile_list = None
for job in jobs[1:]:
    job_id, job_pdb_list, job_paramsfile, job_resfile_list = job.split()
    if job_id and job_pdb_list and job_paramsfile and job_resfile_list:
        job_cmd = shlex.split('qsub python2 '\
                              +fixbb_script+' '\
                              +job_pdb_list+' '\
                              +job_paramsfile+' '\
                              +job_resfile_list)
        print job_cmd
        subprocess.Popen(job_cmd)
    else:
        print 'job sumbission failed at job %s!' % str(job_id)
        sys.exit()
