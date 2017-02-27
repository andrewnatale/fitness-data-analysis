import numpy as np
from Bio import SeqIO
import os
import matplotlib.pyplot as plt
from sklearn import manifold

class distance_calculator(object):

    def __init__(self, filename='BLOSUM62.txt'):
        # load BLOSUM62 matrix as a numpy array
        with open(filename, 'r') as infile:
            lines = infile.readlines()
        # break down the text file into easily indexable array
        pre_matrix = []
        self.indexer = {}
        for line in lines:
            # skip comment lines
            if line[0] == '#':
                continue
            elements = line.split()
            if len(elements) == 24:
                for i,elem in enumerate(elements[0:20]):
                    self.indexer[elem] = i
            if len(elements) == 25:
                if not elements[0] in ['B', 'Z', 'X', '*']:
                    pre_matrix.append([int(n) for n in elements[1:21]])
        # normalize so values range from 0 to 1
        self.matrix = np.array(pre_matrix)
        self.matrix = self.matrix.astype(dtype=np.float64)
        self.matrix -= np.amin(self.matrix)
        self.matrix /= np.amax(self.matrix)
        # invert so that low similarity is large distance
        self.matrix -= 1.0
        self.matrix *= -1.0

    # def lookup(self, aa1, aa2):
    #     aa1 = aa1.upper()
    #     aa2 = aa2.upper()
    #     idx1 = self.indexer[aa1]
    #     idx2 = self.indexer[aa2]
    #     distance = self.matrix[idx1,idx2]
    #     return distance

    def compare_sequences(self, seq1, seq2):
        # check that sequences have the same len
        if len(seq1) != len(seq2):
            print("Can't compare sequences of different lengths!")
            exit(1)
        sum_distance = 0
        for n in range(len(seq1)):
            idx1 = self.indexer[seq1[n].upper()]
            idx2 = self.indexer[seq2[n].upper()]
            sum_distance += self.matrix[idx1,idx2]
        return sum_distance / len(seq1)


class state_tracker(object):

    def __init__(self):
        self.macrostate_count = 0
        self.microstate_count = 0
        self.macrostates = []

    def add_state(self, state_name, filename):
        self.macrostates.append(macrostate(state_name,
                                           self.macrostate_count,
                                           filename,
                                           starting_point=self.microstate_count))
        self.microstate_count += self.macrostates[self.macrostate_count].state_size
        self.macrostate_count += 1

    def report(self):
        print('total macrostates: ', self.macrostate_count)
        print('total microstates: ', self.microstate_count)
        for state in self.macrostates:
            print(state.name, state.id, state.state_size, state.state_range)

class macrostate(object):

    def __init__(self, state_name, state_id, filename, starting_point=0):
        self.name = state_name
        self.id = state_id
        seq_index = starting_point
        self.sequence_dict = {}
        self.sequences = list(SeqIO.parse(filename, 'fasta'))
        for seq in self.sequences:
            self.sequence_dict[seq_index] = str(seq.seq)
            seq_index += 1
        self.state_size = len(self.sequence_dict)
        self.state_range = (starting_point, seq_index)

data_dir = '/Users/anatale/school/UCSF/Kortemme_lab/code/multi-state-design/opt_rnd3'

dc = distance_calculator(filename='BLOSUM62.txt')

tracker = state_tracker()

which_set = 'all'
#which_set = '101_109'
#which_set = '122_127'

tracker.add_state('1i2m', os.path.join(data_dir, '1i2m_%s.fasta' % which_set))
tracker.add_state('1a2k', os.path.join(data_dir, '1a2k_%s.fasta' % which_set))
tracker.add_state('1k5d', os.path.join(data_dir, '1k5d_%s.fasta' % which_set))
tracker.add_state('3gj0', os.path.join(data_dir, '3gj0_%s.fasta' % which_set))
tracker.add_state('importin', os.path.join(data_dir, 'importin_%s.fasta' % which_set))
tracker.add_state('composite', os.path.join(data_dir, 'composite_%s.fasta' % which_set))

tracker.report()

distance_array = np.zeros((tracker.microstate_count, tracker.microstate_count), dtype=np.float64)
# for each sequence...
for stateA in tracker.macrostates:
    for idxA in stateA.sequence_dict:
        sequenceA = tracker.macrostates[stateA.id].sequence_dict[idxA]
        # iterate over every other sequence:
        for stateB in tracker.macrostates:
            for idxB in stateB.sequence_dict:
                sequenceB = tracker.macrostates[stateB.id].sequence_dict[idxB]
                distance_array[idxA,idxB] = dc.compare_sequences(sequenceA, sequenceB)

mds = manifold.MDS(n_components=2, dissimilarity='precomputed', random_state=8)
results = mds.fit(distance_array)
coords = results.embedding_

colors = ['red', 'blue', 'yellow', 'green', 'cyan', 'magenta']

for state in tracker.macrostates:
    print(state.name, colors[state.id])
    plt.scatter(coords[state.state_range[0]:state.state_range[1],0], coords[state.state_range[0]:state.state_range[1],1], color=colors[state.id], marker = 'o')
plt.show()
