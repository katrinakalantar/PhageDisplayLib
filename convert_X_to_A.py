import argparse
from collections import Counter

### create parser object ###
parser = argparse.ArgumentParser(description = 'Generate agilent order')
parser.add_argument('-i','-infile', help='input filepath',type=str)
parser.add_argument('-o','-outfile', help='output filepath',type=str)
args = parser.parse_args()


### utility functions ###
def x_to_alanine(seq):
    if 'X' in seq:
        seq = seq.replace('X','A')
    return(seq)

### initialize variables ###
header = ''
seq = ''

outfile = open(args.o + '.fasta', 'w')

### main ###

with open(args.i,'r') as f:  #open the input file

    for line in f:

        if '>' in line:
            if len(seq) > 0:
                seq = x_to_alanine(seq)
                outfile.write(header + '\n')
                outfile.write(seq + '\n')

            seq = ''     #reset sequence
            header = line.strip()

        else:     #this is not a header line
            seq += line.strip()

    #finish the remaining sequence after last header
    if len(seq) > 0:
        seq = x_to_alanine(seq)
        outfile.write(header + '\n')
        outfile.write(seq + '\n')


outfile.close()