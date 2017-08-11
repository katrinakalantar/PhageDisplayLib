# PhageDisplayLib

### Agilent_Order.py

Adapted from bodonovan's original script.
Generate a .fasta and .txt sequence file containing the agilent order for peptidome.

```
python Agilent_Order.py --help
```


Basic usage:

To generate the agilent order for 62-AA sequences with 31-AA window and STREP and FLAG tags on 5' and 3' ends respectively.

```
python Agilent_Order.py -i input_file.fasta -o output_filename -k 62 -w 31 -fp AGCCATCCGCAGTTCGAGAAA -tp GACTACAAGGACGACGATGAT
```

Special Options:

Case 1: Take in full-length protein sequences and return chopped sequences, still peptide-level AA.

```
python Agilent_Order.py -i input_file.fasta -o output_filename -k 62 -w 31 -p=True
```


Case 2: Take in nucleotide-level sequences (potentially as output from pcr scrambling code) and only add 5' / 3' tags and remove cutsites.

```
pyton Agilent_Order.py -i input_file.fasta -o output_filename -fp five_prime_seq -tp five_prime_seq -r=True
```


### lzw_compress.py

Dependencies:
pip install lzstring

Compress each line in .fasta file to obtain a compression ratio. If ratio < threshold, filter out the sequence.
Returns a .fasta file with low complexity sequences removed and a separate *_removed.fasta file containing all low complexity sequences that failed to pass threshold.

```
python lzw_compress.py -i input_file.fasta -o output_filename -t .5
```

Note: 

