import numpy as np
from collections import Counter
from Bio import SeqIO
from statsmodels import robust
import datetime

class EMPIRICProcessor(object):
    """A class to load and process EMPIRIC style saturation mutagenesis data."""

    # standard aa to # map
    aminotonumber = {'*': 0, 'W': 1,
                     'F': 2, 'Y': 3,
                     'L': 4, 'I': 5,
                     'M': 6, 'V': 7,
                     'C': 8, 'A': 9,
                     'G': 10, 'P': 11,
                     'S': 12, 'T': 13,
                     'N': 14, 'Q': 15,
                     'H': 16, 'R': 17,
                     'K': 18, 'D': 19,
                     'E': 20}

    cutoff = (50,-50)

    def __init__(self):
        """Set triggers to make sure things are done in the right order."""

        self.seq_loaded = False
        self.data_loaded = False
        self.rescaled = False
        self.frequency = False
        # generate a dictionary for reverse number to amino acid lookups
        self.numbertoamino = {}
        for key in self.aminotonumber:
            self.numbertoamino[self.aminotonumber[key]] = key

    def seq_from_fasta(self, native_seq_fasta):
        """Load a reference sequence from a fasta file using Biopython.
        Only uses the first entry in the fasta file if multiple are present."""

        if self.seq_loaded:
            print('Already loaded a reference sequence, cannot load another!')
            exit(0)
        reference = list(SeqIO.parse(native_seq_fasta, 'fasta'))
        self.native_seq = str(reference[0].seq)
        self.header = str(reference[0].id)
        self.seq_loaded = True

    def seq_from_string(self, native_seq_string, identifier=None):
        """Load a reference sequence from an input string."""

        if self.seq_loaded:
            print('Already loaded a reference sequence, cannot load another!')
            exit(0)
        self.native_seq = native_seq_string
        if identifier:
            self.header = identifier
        else:
            self.header = 'null'
        self.seq_loaded = True

    def load_data_csv(self, source_csv):
        """Load preprocessed EMPIRIC data.
        Takes a  a CSV file containing fitness scores and converts it to a numpy
        array. Input file should have a 1 line header and subsequent lines like:
        pos, aa, rep1, rep2, ..., repN
        where 'pos' is an int, 'aa' is a str, and 'reps' are floats."""

        if self.data_loaded:
            print('Already loaded a data file, cannot load another!')
            exit(0)
        with open(source_csv, 'r') as csvfile:
            datalines = csvfile.readlines()
        datalines = [i.strip('\r\n') for i in datalines]
        # remove header of column labels
        datalines.pop(0)
        # read data into a dictionary, replacing the aa letter code with a number
        # which will be used to arrange them (see 'aminotonumber' dictionary)
        fitness_dict = {}
        for line in datalines:
            elem = line.split(',')
            fitness_dict[(elem[0], self.aminotonumber[elem[1]])] \
              = np.array([float(i) for i in elem[2:]])
        # get length of sequence and first resi position
        self.data_length = len(Counter([i[0] for i in fitness_dict]))
        self.first_resi = int(sorted([i[0] for i in fitness_dict])[0])
        # create an empty 3d array of the correct size for all the data
        self.fitness_array = np.zeros((21,self.data_length))
        # fill the array using the data dictionary
        for elem in fitness_dict:
            self.fitness_array[int(elem[1]),int(elem[0])-self.first_resi] \
              = np.mean(fitness_dict[elem])
        # look for bad data points outside of cutoff range
        for idx, j in np.ndenumerate(self.fitness_array):
            if j > self.cutoff[0] or j < self.cutoff[1]:
                self.fitness_array[idx] = np.nan
        self.data_loaded = True

    def generate_labels(self):
        """Compose data labels from available information."""

        if not (self.data_loaded and self.seq_loaded):
            print('Must load sequence and data files before generating labels!')
            exit(0)
        self.sequence_labels = []
        for n in range(self.first_resi, self.first_resi + self.data_length, 1):
            self.sequence_labels.append(('%s%d') % (self.native_seq[n-1], n))
        self.mutation_labels = []
        for key in sorted(self.numbertoamino):
            self.mutation_labels.append(self.numbertoamino[key])

    def rescale(self):
        """Rescale fitness scores so that the max stop codon score becomes zero
        and the average WT amino acid score becomes one."""

        if not self.data_loaded:
            print('No data found! Cannot rescale data before a file has been loaded!')
            exit(0)
        # get max stop codon values
        max_stop = np.amax(self.fitness_array[self.aminotonumber['*'],:])
        # shift everything by the max stop codon value
        self.rescaled_array = self.fitness_array - max_stop
        # set values below zero to zero
        for idx, j in np.ndenumerate(self.rescaled_array):
            if j < 0:
                self.rescaled_array[idx] = 0.0
        # get the average wt residue fitness in the shifted data
        wt_fitness = np.empty(0)
        for j in range(self.first_resi, self.first_resi + self.data_length, 1):
            wt_fitness = np.append(wt_fitness, self.rescaled_array[self.aminotonumber[self.native_seq[j-1]], j - self.first_resi])
        wt_avg = np.mean(wt_fitness[~np.isnan(wt_fitness)])
        # divide everything by the average wt residue fitness
        self.rescaled_array = self.rescaled_array / wt_avg
        self.rescaled = True

    def generate_freqencies(self):
        """Use rescaled fitness scores to calculate a frequency for each
        mutation at each position."""

        if not self.rescaled:
            print('Cannot generate frequencies until data has been rescaled!')
            exit(0)
        self.frequency_array = np.copy(self.rescaled_array)
        for i in range(self.data_length):
            self.frequency_array[:,i] = self.frequency_array[:,i] / np.sum(self.frequency_array[:,i][~np.isnan(self.frequency_array[:,i])])
        self.frequency = True

    def fasta_out_rescaled(self, scale_factor=1000, position_range=None, outfile='EMPIRIC_to_seq.fasta'):
        """Write out frequency data for a range of positions as a bunch of
        dummy fasta format sequences which can be used to reconstruct the
        frequency data (with some loss of precision)."""

        if not self.frequency:
            print('Cannot generate fasta file before frequency calculation!')
            exit(0)
        if not position_range:
            selected_positions = range(self.first_resi, self.first_resi + self.data_length, 1)
        else:
            selected_positions = [int(i) for i in position_range]
        # empty char array to hold sequences as they are constructed
        seq_array = np.zeros((len(selected_positions), scale_factor), dtype=str)
        for pidx,position in enumerate(selected_positions):
            # get array index
            for lidx,label in enumerate(self.sequence_labels):
                if int(label[1:]) == position:
                    target_idx = lidx
            # first scale up the floats and convert NaN to zero
            unrounded = np.nan_to_num(self.frequency_array[1:,target_idx]) * scale_factor
            # round and convert to ints
            rounded = np.around(unrounded).astype(int)
            # adjust highest value so that everything sums correctly
            diff = np.sum(rounded) - scale_factor
            if diff != 0:
                rounded[np.argmax(rounded)] -= diff
            # write frequecies into array as amino acids
            seq_idx = 0
            for idx,j in np.ndenumerate(rounded):
                aa_count = 0
                while aa_count < j:
                    seq_array[pidx,seq_idx] = str(self.numbertoamino[idx[0]+1])
                    aa_count += 1
                    seq_idx += 1
        # write from array into output file
        with open(outfile, 'w') as fasta_out:
            for i in range(scale_factor):
                fasta_out.write('>seq %d\n' % (i+1))
                fasta_out.write('%s\n' % ''.join(seq_array[:,i]))

    def errorstats(self):
        """Use the spread in WT fitness scores to estimate an error on each score"""

        if not self.data_loaded:
            print('No data found! Cannot calculate error before a file has been loaded!')
            exit(0)
        # using only unscaled data:
        self.stop_vals = self.fitness_array[self.aminotonumber['*'],:]
        self.wt_vals = np.empty(0)
        for j in range(self.first_resi, self.first_resi + self.data_length, 1):
            self.wt_vals = np.append(self.wt_vals, self.fitness_array[self.aminotonumber[self.native_seq[j-1]], j - self.first_resi])
        self.stop_vals = self.stop_vals[~np.isnan(self.stop_vals)]
        self.wt_vals = self.wt_vals[~np.isnan(self.wt_vals)]
        self.stop_mean = np.mean(self.stop_vals)
        self.stop_std = np.std(self.stop_vals)
        self.stop_mad = robust.scale.mad(self.stop_vals)
        self.wt_mean = np.mean(self.wt_vals)
        self.wt_std = np.std(self.wt_vals)
        self.wt_mad = robust.scale.mad(self.wt_vals)
        print('stop')
        print('mean %f, std %f, mad %f' % (self.stop_mean, self.stop_std, self.stop_mad))
        print(np.amax(self.stop_vals), np.amin(self.stop_vals))
        print('wt')
        print('mean %f, std %f, mad %f' % (self.wt_mean, self.wt_std, self.wt_mad))
        print(np.amax(self.wt_vals), np.amin(self.wt_vals))
