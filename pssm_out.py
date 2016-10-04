#!/usr/bin/env python2
#
# script to write normalized data into a position specific scoring matrix
# in the transfac file format that can be processed by the weblogo
# application to generate a sequence logo

import numpy as np
from fitness_analysis import process_fitness
from fasta_aa import load_fasta_aa

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
        # write data for each position on a line
        # start by writing the position labels
        position_count = 0
        for i in seq_range:
            pssm_file.write('%s ' % i[1:])
            for j in data[:,position_count]:
                pssm_file.write('%s ' % str(j))
            pssm_file.write('\n')
            position_count += 1
        # write footer
        pssm_file.write('XX\n')
        pssm_file.write('//\n')

# hardcoded data filenames
native_seq_file = '../../logo/P32835.fasta'
fitness_data_file = '../../logo/gsp_log_data_101-140.csv'

a, b, c = process_fitness(fitness_data_file, native_seq_file, 'stop')

pssm_out(a[1:,:], b, c[1:], 'header', 'pssm_stop.txt')
