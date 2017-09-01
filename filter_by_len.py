
import sys
import Bio
from Bio import SeqIO

input_file = open(sys.argv[1],'r').readlines()
output = []
remove = []

header = None
seq = None

for line in input_file:
    if '>' in line:
        header = line.strip()
    else:
        seq = line.strip()

    if header and seq:
        #do stuff

        if len(seq) < 30:
            remove.append(header)
            remove.append(seq)
            header = None
            seq = None
        else:
            output.append(header)
            output.append(seq)
            header = None
            seq = None

 
keep_output = '.'.join(sys.argv[1].split('.')[0:-1]) + '.keep.fasta'
removed_output = '.'.join(sys.argv[1].split('.')[0:-1]) + '.remove.fasta'

output_file = open(keep_output, 'w')
output_file.write('\n'.join(output))

removed_file = open(removed_output, 'w')
removed_file.write('\n'.join(remove))
