# write out a pymol script to load all the pdbs in a directory

import os

with open('load.pml', 'w') as pymol_sc:
    pymol_sc.write('reinit\n')
    for filename in os.listdir(os.getcwd()):
        if filename.endswith('.pdb'):
            pymol_sc.write('load %s\n' % filename)
    pymol_sc.write('hide everything\n')
    pymol_sc.write('show ribbon, chain A\n')
    pymol_sc.write('show sticks, chain X\n')
