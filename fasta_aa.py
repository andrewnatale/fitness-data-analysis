#!/usr/bin/env python2
# load a fasta format amino acid sequence
# into a dictionary keyed by reside number

def load_fasta_aa(filename):
    # load file and strip whitespace
    fasta_file = open(filename, 'r')
    fasta_lines = fasta_file.readlines()
    fasta_file.close()
    fasta_lines = [i.strip('\n') for i in fasta_lines]
    # remove header
    header = fasta_lines.pop(0)
    # join lines into a string
    fasta_string = ''.join(fasta_lines)
    # build a dictionary keyed by sequence number
    fasta_dict = {}
    fasta_dict['header'] = header
    resi_no = 1
    for resi in fasta_string:
        fasta_dict[resi_no] = resi
        resi_no += 1
    return fasta_dict

#print load_fasta('P62825.fasta')
