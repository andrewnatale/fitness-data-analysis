#!/usr/bin/env python2
#
# script to write normalized data into a position specific scoring matrix
# in the transfac file format that can be processed by the weblogo
# application to generate a sequence logo

import numpy as np
from fitness_analysis import process_fitness

def pssm_out(data, seq_range, positions, identifier, out_file):
    # open a file for writing
    with open(out_file, 'w') as pssm_file:
        # put an identifier on the first line
        pssm_file.write('ID %s\n' % identifier)
        # write position values
        pssm_file.write('P0 ')
        for elem in positions:
            pssm_file.write('%s ' % elem)
        pssm_file.write('\n')
        # write data for each position
        #for elem in seq_range:
        #    pass
        # write footer
        pssm_file.write('XX\n')
        pssm_file.write('//\n')

# hardcoded data filenames
native_seq_file = 'P32835.fasta'
fitness_data_file = 'gsp_log_data_101-140.csv'

a, b, c = process_fitness(fitness_data_file, native_seq_file)

pssm_out(None, None, column_labels, native_seq['header'], 'stuff.test')
