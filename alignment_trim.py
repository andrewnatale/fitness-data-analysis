from Bio import SeqIO

infile = 'Ran_uniprot_aln.fasta'

alignment = list(SeqIO.parse(infile, 'fasta'))

print(len(alignment))

trim = False
if trim:
    outfile = 'Ran_aln_emp40.fasta'
    with open(outfile, 'w') as fastafile:
        for seq in alignment:
            fastafile.write('>'+seq.id+'\n')
            #fastafile.write(str(seq.seq[441:444])+str(seq.seq[445:451])+'\n')
            #fastafile.write(str(seq.seq[490:493])+str(seq.seq[498:501])+'\n')
            fastafile.write(str(seq.seq[441:444])+
                            str(seq.seq[445:452])+
                            str(seq.seq[455:459])+
                            str(seq.seq[461:463])+
                            str(seq.seq[474:479])+
                            str(seq.seq[490:493])+
                            str(seq.seq[498:503])+
                            str(seq.seq[520:524])+
                            str(seq.seq[526:533])+
                            '\n')
