import random
import os
import ast

aminotonumber = {'W': 1,
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

numbertoamino = {}
for key in aminotonumber:
    numbertoamino[aminotonumber[key]] = key

def write_random_fasta(outfile, entries=10, seq_length=40):
    n = 1
    m = 1
    with open(outfile, 'w') as fasta:
        while n <= entries:
            fasta.write('>sequence %d\n' % n)
            seq_string = ''
            while m <= seq_length:
                seq_string += numbertoamino[random.randint(1,20)]
                m += 1
            fasta.write('%s\n' % seq_string)
            m = 1
            n +=1

#write_random_fasta('test.fasta')

def write_random_tsv(outfile, mstates=['STATE-A', 'STATE-B', 'STATE-C'], seq_length=40):
    with open(outfile, 'w') as tsvfile:
        tsvfile.write('dummy header\n')
        for state in mstates:
            for n in range(seq_length):
                tsvfile.write('%s\t' % state)
                tsvfile.write('0.9\t')
                tsvfile.write('60\t')
                tsvfile.write('mean\t')
                tsvfile.write('%s\t' % str(n+1))
                tsvfile.write('{')
                for amino in aminotonumber:
                    tsvfile.write("'%s':%s," % (amino, str(random.uniform(-100.,100.))))
                tsvfile.write('}\n')

#write_random_tsv('data_test.tsv')

def test_tsv_parse(infile):
    firstline = True
    with open(infile, 'r') as tsvfile:
        for line in tsvfile:
            if firstline:
                firstline = False
            else:
                entry = line.split('\t')
                macrostate = entry[0];		# string
                backrubT = entry[1];			# string
                ensembleS = entry[2];			# string
                boltzmanT = entry[3];			# string
                position = entry[4];			# string
                energies = ast.literal_eval(entry[5]);
                print(macrostate, backrubT, ensembleS, boltzmanT, position)
                print(energies['W'])

test_tsv_parse('testing.tsv')
