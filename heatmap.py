#!/usr/bin/env python2
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from fitness_analysis import process_fitness

def plot_hmap(data, row_labels, column_labels):
    sns.set(font_scale=0.7)

    grid_kws = {"height_ratios": (.9, .05), "hspace": .3}
    f, (ax, cbar_ax) = plt.subplots(2, gridspec_kw=grid_kws)
    ax = sns.heatmap(data, ax=ax,
                     cbar_ax=cbar_ax,
                     cbar_kws={"orientation": "horizontal"},
                     xticklabels=row_labels,
                     yticklabels=column_labels)
    plt.sca(ax)
    plt.yticks(rotation=0)
    plt.xticks(rotation=90)
    plt.show()

# hardcoded data filenames
native_seq_file = 'P32835.fasta'
fitness_data_file = 'gsp_log_data_101-140.csv'

a, b, c = process_fitness(fitness_data_file, native_seq_file)

# slice off stop codon data points
plot_hmap(a[1:,:], b, c[1:])
