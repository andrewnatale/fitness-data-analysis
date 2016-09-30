#!/usr/bin/env python2
import numpy as np
from fitness_analysis import process_fitness

native_seq_file = '../../logo/P32835.fasta'
fitness_data_file = '../../logo/gsp_log_data_101-140.csv'

a, b, c = process_fitness(fitness_data_file, native_seq_file)
