#!/usr/bin/env python2
# test rosetta features reporting

import os
import sys
import shlex
import subprocess

rosetta_rscripts = '/Users/anatale/Rosetta/main/source/bin/rosetta_scripts.default.macosclangrelease'
rosetta_db_dir = "/Users/anatale/Rosetta/main/database"

structures_path = os.getcwd()
#structures_list_file = ''
xml_script_file = sys.argv[1]
params_file = '1K5D_TP0_backrub.params'

db_out_name = 'test.db3'

listfile_name = 'pdb_list.lst'
with open(listfile_name, 'w') as listfile:
    for filename in os.listdir(structures_path):
        if filename.endswith('.pdb'):
            listfile.write('%s\n' % filename)

rscripts_cmd = shlex.split(rosetta_rscripts+
  ' -database '+rosetta_db_dir+
  ' -l '+listfile_name+\
  ' -extra_res_fa '+params_file+
  ' -out:use_database '
  #' -inout:dbms:database_mode sqlite3'+
  #' -inout:database_filename='+db_out_name+
  ' -parser:protocol '+xml_script_file+
  ' -out:nooutput')

print rscripts_cmd

process = subprocess.Popen(rscripts_cmd)
returncode = process.wait()
