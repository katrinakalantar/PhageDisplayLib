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


### homopolymer_removal.py

Remove sequences with strings of a single amino acid of length K, where default K = 8

```
python homopolyme_removal.py -i input_file.fasta -o output_filename -k 8
```

To determine the appropriate threshold, make use of scripts in /test_scripts. The lzw_compress.TEST.py file (run in the same way as lzw_compress.py) outputs a .csv file containing sequences and headers sorted by compression ratio. This file is referenced as input to the .ipynb jupyter notebook. NOTE: will need to update file path to the appropriate path for the .csv file. The jupyter notebook allows for visualization of the distribution of both homopolymer occurrence and compression ratio.


### convert_X_to_A.py

Loops through all sequences and converts all 'X' amino acides to 'A'.

```
python convert_X_to_A.py -i input_file.fast -o output_filename
```


