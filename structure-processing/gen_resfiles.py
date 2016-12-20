#!/usr/bin/env python2

# Generate resfiles

from Bio.PDB import *
import os
import pandas as pd
import sys

# configurable script parameters:

# target residues as a list of ints
targets_yeast = [101,102,103,105,106,107,108,109,122,123,124,125,126,127]
# radius (angstroms) around the target residues to look for res to repack
radius = 9
# output directory
output_prefix = '/Users/anatale/Documents/school/UCSF/Kortemme_lab/rosetta_staging/resfiles'
index_file = '/Users/anatale/Documents/school/UCSF/Kortemme_lab/code/structure-processing/aux_files/Gsp1_index_all.txt'

amino_acids = 'WFYLIMVCAGPSTNQHRKDE'

# make output directory
try:
    os.makedirs(output_prefix)
except OSError:
    if not os.path.isdir(output_prefix):
        raise

# load the index file into a dataframe
seq_index = pd.read_table(index_file)

# setup biopython tools
parser = PDBParser()
io = PDBIO()

def gen_resfiles(structure, target_number):
    # translate yeast seq# to pdb seq# using index file
    pdb_num = None
    #print structure.id[:-8].lower(), target_number
    for index, row in seq_index.iterrows():
        if row['pdb'] == structure.id[:-8].lower():
            if row['yeast_seq_num'] == target_number:
                pdb_num = row['pdb_res_num']
    #print 'pdb_num:', pdb_num
    chainA = structure[0]['A']
    # use all the atoms in the target residue as search terms
    target_atoms = Selection.unfold_entities(chainA[pdb_num], 'A')
    # define the search area; all atoms except those in the target res and in the ligand
    atom_list = [atom for atom in structure.get_atoms() \
                 if not ((atom.get_parent().id[1] == pdb_num and atom.get_parent().get_parent().id == 'A')\
                 or atom.get_parent().get_parent().id == 'X')]
    ns = NeighborSearch(atom_list)
    # search and create a set of repacking residues
    shell = {(res.id[1], res.get_parent().id) for target_atom in target_atoms \
      for res in ns.search(target_atom.coord, radius, 'R')}
    # write a resfile
    resfile_list = []
    for aa in amino_acids:
        with open('%s_%s%s.res' %(structure.id, str(pdb_num), aa), 'w') as resfile:
            resfile.write('NATRO\n')
            resfile.write('START\n')
            resfile.write('%s A PIKAA %s\n' % (pdb_num, aa))
            for res in sorted(shell):
                resfile.write('%s %s NATAA\n' % (res[0], res[1]))
        resfile_list.append('%s_%s%s.res' %(structure.id, str(pdb_num), aa))
    return resfile_list

# load all pdb files in cwd into a list
pdb_files = []
for pdb_file in os.listdir(os.getcwd()):
    if pdb_file.endswith(".pdb"):
        pdb_files.append(parser.get_structure(pdb_file[:-4], pdb_file))

# we expect pdb files with '_backrub' appended, while the index file does not
# so slice it off to check if the file is indexed
for obj in pdb_files:
    print obj.id
    if obj.id.lower()[:-8] in seq_index['pdb'].tolist():
        print 'indexed'
    else:
        print 'not indexed'

# make resfiles for each pdb in its own dir
for structure in pdb_files:
    # make a directory for the resfiles
    print 'generating resfiles: ', structure.id
    try:
        os.makedirs(os.path.join(output_prefix, structure.id))
    except OSError:
        if not os.path.isdir(os.path.join(output_prefix, structure.id)):
            raise
    os.chdir(os.path.join(output_prefix, structure.id))
    resfiles_for_structure = []
    for target in targets_yeast:
        resfiles_for_target = gen_resfiles(structure, target)
        resfiles_for_structure.extend(resfiles_for_target)
    with open('%s_resfiles.lst' % structure.id, 'w') as resfile_list_file:
        resfile_num = 1
        for resfile in resfiles_for_structure:
            resfile_list_file.write('%s %s\n' % (str(resfile_num), resfile))
            resfile_num += 1
    os.chdir(output_prefix)
