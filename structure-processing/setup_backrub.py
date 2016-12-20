#!/usr/bin/env python2

# find pivot residues and write params files for GNP ligand containing structures

# !!!make sure to sanity check all the output files -
# this script may not handle bad input gracefully!!!

from Bio.PDB import *
import re
import subprocess
import shlex
import os
import pandas as pd
import sys

# configurable script parameters:

# target residues for defining the pivot residues, as a set of ints
targets_yeast = {122,123,124,125,126,127}
# radius (angstroms) around the target residues to look for pivot residues
radius = 9
# output directory
output_dir = 'backrub_ready'
output_prefix = os.path.join(os.getcwd(), output_dir)
output_suffix = 'backrub'
#output_prefix = '/Users/anatale/Documents/school/UCSF/Kortemme_lab/code/structure_scripts/pdbs/test'
# needs this script to actually write the params file
rosetta_pyscripts = ''
molfile_to_params = 'molfile_to_params.py'
# generate pml script?
pml_file = True

# make output directory
try:
    os.makedirs(output_prefix)
except OSError:
    if not os.path.isdir(output_prefix):
        raise
# load the index file into a dataframe
seq_index = pd.read_table(sys.argv[1])

# setup biopython tools
parser = PDBParser()
io = PDBIO()

find_chi = re.compile('CHI')

# find ligands, the Bio.PDB way
find_gnp = re.compile('H_G..')
class GNPSelect(Select):
    def accept_residue(self, residue):
        if find_gnp.match(residue.id[0]):
            return 1
        else:
            return 0

# clean structures - remove water and unneeded ligands (like SO4)
# based on this list - update as needed
reject_res = ['W', 'H_SO4']
class NoH2OSelect(Select):
    def accept_residue(self, residue):
        if residue.id[0] in reject_res:
            return 0
        else:
            return 1

class PivotsParamsGen(object):

    def __init__(self, pdb_file):
        self.structure = parser.get_structure(pdb_file[:-4], pdb_file)
        # this must be lowercase to properly address the index file
        self.name = self.structure.id.lower()
        if self.name in seq_index['pdb'].tolist():
            self.indexed = True
        else:
            self.indexed = False

    def generate_pivots(self):
        self.targets_pdb_num = set()
        # translate yeast seq# to pdb seq# using index file
        for index, row in seq_index.iterrows():
            if row['pdb'] == self.name:
                if row['yeast_seq_num'] in targets_yeast:
                    self.targets_pdb_num.add(row['pdb_res_num'])
        # Gsp1/RAN is always chain A
        self.chainA = self.structure[0]['A']
        # use all the atoms in the target residues as search terms
        self.targets = [self.chainA[res] for res in self.targets_pdb_num]
        self.target_atoms = Selection.unfold_entities(self.targets, 'A')
        # define the search area; for now all backbone and CB heavy atoms
        self.atom_list = [atom for atom in self.chainA.get_atoms() \
                          if atom.name == 'CA' \
                          or atom.name == 'C' \
                          or atom.name == 'N' \
                          or atom.name == 'CB']
        self.ns = NeighborSearch(self.atom_list)
        # search and create a set of pivot residues
        self.pivots_pdb = {res.id[1] for target_atom in self.target_atoms \
          for res in self.ns.search(target_atom.coord, radius, 'R')}
        # convert pdb numbering to yeast numbering
        self.pivots_yeast = set()
        for index, row in seq_index.iterrows():
            if row['pdb'] == self.name:
                if row['pdb_res_num'] in self.pivots_pdb:
                    self.pivots_yeast.add(row['yeast_seq_num'])
        #print sorted(self.targets_pdb_num)
        #print sorted(self.pivots_pdb)
        #print sorted(self.pivots_yeast)

    # take a set of pivot res in yeast# and return rosetta#
    def retranslate(self, pivot_set):
        # translate yeast pivots to pdb pivots
        self.pivots_pdb_retr = set()
        for index, row in seq_index.iterrows():
            if row['pdb'] == self.name:
                if row['yeast_seq_num'] in pivot_set:
                    self.pivots_pdb_retr.add(row['pdb_res_num'])
        # get a list of res# in pdb file
        self.pdb_res_list = [res.id[1] for res in self.chainA.get_residues() \
                             if res.id[0] == ' ']
        # make a list of rosetta pivots
        self.pivots_rosetta = [(index+1) for index, res in enumerate(self.pdb_res_list) \
                               if res in self.pivots_pdb_retr]

    def generate_params(self):
        # identify gnp ligand
        self.gnp_res_name = None
        for res in Selection.unfold_entities(self.structure, 'R'):
            if find_gnp.match(res.id[0]):
                self.gnp_res_name = res.id[0][2:]
        # save a cleaned up structure and a separate pdb with only ligand
        io.set_structure(self.structure)
        io.save(os.path.join(output_prefix, '%s_%s.pdb' % (self.structure.id, output_suffix)), NoH2OSelect())
        # if there is no apparent ligand stop now
        if not self.gnp_res_name:
            print '****************************'
            print 'No ligand found for this PDB'
            print '****************************'
            self.new_params_filename = 'NA'
        # if there is a ligand, convert it to a .params file
        else:
            self.lig_pdb_file = '%s_GNP.pdb' % self.structure.id
            io.save(self.lig_pdb_file, GNPSelect())
            # use babel to convert stripped pdb to mol
            self.lig_mol_file = '%s.mol' % self.structure.id
            subprocess.call(shlex.split('babel -ipdb '+self.lig_pdb_file+\
              ' -omol '+self.lig_mol_file+' -h --title '+self.gnp_res_name))
            # now use a rosetta script to convert that mol file to a params file
            subprocess.call(shlex.split('python2 '\
              +os.path.join(rosetta_pyscripts, molfile_to_params)+\
              ' '+self.lig_mol_file+' --name '+self.gnp_res_name))
            # clean and rename params file
            with open('%s.params' % self.gnp_res_name, 'r') as params:
                self.params_lines = params.readlines()
            self.new_params_filename = '%s_%s.params' % (self.structure.id, output_suffix)
            # find chi angle descriptors in the params file and delete
            with open(os.path.join(output_prefix, self.new_params_filename), 'w') as params:
                for line in self.params_lines:
                    if not find_chi.search(line):
                        params.write(line)
            # clean up junk files
            os.remove('%s.params' % self.gnp_res_name)
            os.remove('%s_0001.pdb' % self.gnp_res_name)
            os.remove('%s_GNP.pdb' % self.structure.id)
            os.remove(self.lig_mol_file)

# load all pdb files in cwd into PivotsParamsGen objects
pdb_files = []
for filename in os.listdir(os.getcwd()):
    if filename.endswith(".pdb"):
        pdb_files.append(PivotsParamsGen(filename))

for obj in pdb_files:
    print obj.indexed

# generate an inclusive set of pivots from the sets for each structure
yeast_pivots_all = set()
for obj in pdb_files:
    if not obj.indexed:
        continue
    obj.generate_pivots()
    print sorted(obj.pivots_yeast)
    yeast_pivots_all.update(obj.pivots_yeast)
print sorted(yeast_pivots_all)

# fill in gaps of up to 2 res in the set
for res in sorted(yeast_pivots_all):
    if res + 2 in yeast_pivots_all:
        yeast_pivots_all.add(res + 1)
    elif res + 3 in yeast_pivots_all:
        yeast_pivots_all.add(res + 1)
        yeast_pivots_all.add(res + 2)
    if res - 2 in yeast_pivots_all:
        yeast_pivots_all.add(res - 1)
    elif res - 3 in yeast_pivots_all:
        yeast_pivots_all.add(res - 1)
        yeast_pivots_all.add(res - 2)

# filter the set so that it only contains contiguous regions of 2 or more res
for res in sorted(yeast_pivots_all):
    if not (res + 1 in yeast_pivots_all or res - 1 in yeast_pivots_all):
        yeast_pivots_all.remove(res)

for obj in pdb_files:
    if not obj.indexed:
        continue
    obj.retranslate(yeast_pivots_all)

if pml_file:
    with open(os.path.join(output_prefix, 'check_pivots.pml'), 'w') as pml:
        pml.write('reinit\n')
        for obj in pdb_files:
            if not obj.indexed:
                continue
            pml.write('load %s_%s.pdb\n' % (obj.structure.id, output_suffix))
            pml.write('hide everything, %s\n' % obj.structure.id)
            pml.write('show car, chain A\n')
            for pivot_res in sorted(obj.pivots_pdb_retr):
                pml.write('color red, %s and chain A and resi %d\n' % \
                  (obj.structure.id, pivot_res))
            for target_res in sorted(obj.targets_pdb_num):
                pml.write('show sticks, %s and chain A and resi %d\n' % \
                  (obj.structure.id, target_res))

# generate a .params file and find pivot residues for every pdb file
# in the cwd, then write a file called params_list.txt with the info
with open(os.path.join(output_prefix, 'params_list.txt'), 'a+') as params_list:
     # write a header on line 1
     params_list.write('pdb_file params_file pivot_residues\n')
     # generate params files and pivot residue info,
     # then write it in a format the rosetta job submission script can read
     for obj in pdb_files:
         if not obj.indexed:
             continue
         obj.generate_params()
         params_list.write('%s_%s.pdb %s %s\n' % \
           (obj.structure.id, output_suffix, obj.new_params_filename, \
           ' '.join([str(i) for i in obj.pivots_rosetta])))
