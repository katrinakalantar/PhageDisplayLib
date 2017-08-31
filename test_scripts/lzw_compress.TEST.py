import argparse
import lzstring
import pandas as pd

### create parser object ###
parser = argparse.ArgumentParser(description = 'Generate agilent order')
parser.add_argument('-i','-infile', help='input filepath',type=str)
parser.add_argument('-o','-outfile', help='output filepath',type=str)
parser.add_argument('-t','-threshold', help='compression ratio threshold; default = .5', default=.5, nargs='?', type=float)
args = parser.parse_args()


### utility functions ###
def utf8len(s):
    return len(s.encode('utf-8'))

def get_compression_ratio(seq):
    x = lzstring.LZString()
    compressed = x.compressToBase64(seq)
    ratio = utf8len(compressed) / utf8len(seq)

    #if ratio < .6:
    #    print(seq)
    #    print(ratio)

    return ratio

### initialize variables ###
header = ''
seq = ''

outfile = open(args.o + '.fasta', 'w')
outfile_removed = open(args.o + '_removed.fasta', 'w')


dictionary_of_seqs = {'seq':[],'id':[],'ratio':[]}


### main ###

with open(args.i,'r') as f:  #open the input file

    for line in f:

        if '>' in line:
            if len(seq) > 0:
                cr = get_compression_ratio(seq)
                if cr > args.t:  #if compression ratio > the specified threshold
                    outfile.write(header + '\n')
                    outfile.write(seq + '\n')

                    dictionary_of_seqs['seq'].append(seq)
                    dictionary_of_seqs['id'].append(header)
                    dictionary_of_seqs['ratio'].append(cr)

                else:
                    outfile_removed.write(header + '\n')
                    outfile_removed.write(seq + '\n')

                    dictionary_of_seqs['seq'].append(seq)
                    dictionary_of_seqs['id'].append(header)
                    dictionary_of_seqs['ratio'].append(cr)

            seq = ''     #reset sequence
            header = line.strip()

        else:     #this is not a header line
            seq += line.strip()

    #finish the remaining sequence after last header
    if len(seq) > 0:
        cr = get_compression_ratio(seq)
        if cr > args.t:  #if compression ratio > the specified threshold
            outfile.write(header + '\n')
            outfile.write(seq + '\n')


df = pd.DataFrame.from_dict(dictionary_of_seqs)
df.sort_values('ratio', inplace=True)
df = df[['ratio','seq','id']]
df.to_csv('output_dataframe.csv')