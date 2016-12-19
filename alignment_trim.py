from Bio import SeqIO

infile = 'Ran_uniprot_aln.fasta'

alignment = list(SeqIO.parse(infile, 'fasta'))

outfile = 'Ran_aln_trim.fasta'
with open(outfile, 'w') as fastafile:
    for seq in alignment:
        fastafile.write('>'+seq.id+'\n')
        fastafile.write(str(seq.seq[490:493])+str(seq.seq[498:501])+'\n')
