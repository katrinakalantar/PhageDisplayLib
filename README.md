#
PhageDisplayLib

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


