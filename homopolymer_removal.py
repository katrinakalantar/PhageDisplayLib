import argparse
from collections import Counter

### create parser object ###
parser = argparse.ArgumentParser(description = 'Generate agilent order')
parser.add_argument('-i','-infile', help='input filepath',type=str)
parser.add_argument('-o','-outfile', help='output filepath',type=str)
parser.add_argument('-k','-kmer_len', help='kmer length to remove; default = 8', default=8, nargs='?', type=int)
args = parser.parse_args()


### utility functions ###
def check_for_homopolymer(seq,k,verbose=False):
    is_homopolymer = False
    kmer_dict = Counter([seq[i:i+k] for i in range(len(seq)-k)])
    for k in kmer_dict.keys():
        if len(set(k)) == 1:
            is_homopolymer = True  
            
    if verbose:
        if is_homopolymer:
            print(seq)
            
        if not is_homopolymer:
            if 'X' in seq:
                print('contains Xs but not homopolymer: ' + seq)
    return(is_homopolymer)

### initialize variables ###
header = ''
seq = ''

outfile = open(args.o + '.fasta', 'w')
outfile_removed = open(args.o + '_removed.fasta', 'w')

### main ###

with open(args.i,'r') as f:  #open the input file

    for line in f:

        if '>' in line:
            if len(seq) > 0:
                is_homopolymer = check_for_homopolymer(seq,args.k)
                if is_homopolymer:  #if sequence is a homopolymer, save to "removed" file
                    outfile_removed.write(header + '\n')
                    outfile_removed.write(seq + '\n')
                else:
                    outfile.write(header + '\n')
                    outfile.write(seq + '\n')
                    

            seq = ''     #reset sequence
            header = line.strip()

        else:     #this is not a header line
            seq += line.strip()

    #finish the remaining sequence after last header
    if len(seq) > 0:
        is_homopolymer = check_for_homopolymer(seq,args.k)
        if is_homopolymer:  #if sequence is a homopolymer, save to "removed" file
            outfile_removed.write(header + '\n')
            outfile_removed.write(seq + '\n')
        else:
            outfile.write(header + '\n')
            outfile.write(seq + '\n')


outfile.close()
outfile_removed.close()