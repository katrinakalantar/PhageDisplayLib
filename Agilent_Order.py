
# coding: utf-8

# # Agilent order for display library
# Agilent MTA came through and the prices are pretty decent. 170 and 190mers are the options. Need to take a look at the optimal lengths and compression parameters to get the best and most cost effective library.
#
# Substantial modifications made by KKalantar on 3/27

# In[63]:

import pandas as pd
import json
import seaborn as sns
import matplotlib as plt
import Bio
from Bio.Seq import Seq
from collections import defaultdict
import json
import subprocess
import sys
import os
import warnings
import Bio
from Bio.Seq import Seq
import argparse
import math



###Argument Parsing Functions###

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

parser = argparse.ArgumentParser(description = 'Generate agilent order')
parser.add_argument('-i','-infile', help='input filepath',type=str)
parser.add_argument('-o','-outfile', help='output filepath (root name, no extension)',type=str)
parser.add_argument('-k','-kmer_size', help='length of kmer; default = 62', default=62, nargs='?', type=int)
parser.add_argument('-w','-window_size', help='length of window; default = 31', default=31, nargs='?', type=int)
parser.add_argument('-fp','-fiveprime', help='five prime end barcode sequence (nucleotide); default None',nargs='?',default='', type=str)  #, const='AGCCATCCGCAGTTCGAGAAA')  #default = STREP
parser.add_argument('-tp','-threeprime', help='three prime end barcode sequence (nucleotide); default None',nargs='?',default='', type=str)  #, const='GACTACAAGGACGACGATGAT')
parser.add_argument('-p','-convert2peptide', help='return chopped up peptide sequences', type=str2bool, nargs='?',const=False)
parser.add_argument('-r','-restrictiononly', type=str2bool, nargs='?', const=False,
help="Run only restriction checking, default = False \n *note: -r=True is required when working with NT input")
args = parser.parse_args()
print(args)


# In[81]:

###Functions###


def strip_stop_codon(sequence):
    '''strips stop codon and \n character at the end of each sequence'''
    sequence = sequence.strip('\n')
    if sequence[-3:]  == 'TAA'  or sequence[-3:] == 'TAG' or sequence[-3:] == 'TGA' :
    	sequence = sequence[0:-3]
    return sequence


def split_kmers(sequence, kmer_size, window_size):
    '''given a sequence, this will return a list kmers of a given length based on a given sliding window'''
    ##Handles end/edge cases where the length of the sequence isn't evenly divisible by the sliding window size
    kmer_list = []
    short_list = []
    if len(sequence) < kmer_size:
        #short_list.append(sequence)
        #print(sequence)
        pad_len = kmer_size - len(sequence)
        if(pad_len%2==1):
            pad1 = int(math.floor(pad_len/2))
            pad2 = int(math.ceil(pad_len/2))
        else:
            pad1 = int(pad_len/2)
            pad2 = int(pad_len/2)
        sequence = ''.join(['A']*pad1) + sequence + ''.join(['A']*pad2)
    split_seqs = [sequence[i:i+kmer_size] for i in range(0, len(sequence), window_size)]
    split_seqs.append(sequence[-kmer_size:]) # Ensure that you get the last kmer
    for kmer in split_seqs:
        if len(kmer) == kmer_size and kmer not in kmer_list:
            kmer_list.append(kmer)

    return kmer_list


def split_kmer_dict(sequence_dictionary, kmer_size, window_size):
    '''for every sequence in dictionary, split into kmers of size kmer_size with overlap window_size'''

    split_kmer_dict = {}
    for key in sequence_dictionary.keys():
        sequences = split_kmers(sequence_dictionary[key], kmer_size, window_size)
        split_kmer_dict[key] = sequences

    return(split_kmer_dict)

def open_infile(infile):
    '''opens the fasta file, reads in  peptide sequences, reverse-tranlates protein --> DNA, returns list of DNA sequences'''

    peptide_dict = {}

    with open(infile, 'r') as _in:
        data = _in.readlines()
        for line in data:
            if line.startswith('>'):
                seq_ID = line.strip('\n')
                peptide_dict[seq_ID] = ''
            else:
                peptide = line.strip('\n')
                peptide_dict[seq_ID]+=peptide

    return(peptide_dict)


def prot_to_dna(protein):
    '''Translates protein sequence to e. coli preferred codons'''

    E_coli_codons = {'X': 'GCT', 'U': 'TGC', 'F':'TTT', 'L':'CTG', 'I': 'ATT', 'M':'ATG', 'V':'GTG', 'S':'TCT', 'P':'CCG', 'T':'ACC', 'A':'GCG', 'Y':'TAT', '*':'TAA', 'H':'CAT', 'Q':'CAG', 'N': 'AAC' ,'K':'AAA', 'D':'GAT', 'E':'GAA', 'C':'TGC', 'W':'TGG', 'R':'CGT' , 'S':'AGC', 'G':'GGC' , 'B':'GAT' , 'Z':'GAA' , 'J':'CTG'}
    #handles selenocystein and X's with alanines... no other ns-AAs handled
    DNA = ''
    for amino in protein.strip('\n'):
        try:
            DNA = DNA + E_coli_codons[amino]
        except:
#             unnatural_AAs.append(amino)
            print ('unnatural amino!?' + ' ' + amino)
    return DNA


def generate_DNA_dict(peptide_dict):
    '''pass me a dictionary of seq_IDs:peptides and I will give you a dictionary of seq_IDs:e.coli optimized oligos'''

    oligo_dict = {}
    for seq_ID in peptide_dict.keys():
        oligo_dict[seq_ID] = prot_to_dna(peptide_dict[seq_ID])

    return oligo_dict


def sanity_check(oligo_dict, peptide_dict):
    '''checks if the reverse_translated sequences in the oligo dict do, in fact, code for the peptide'''

    for seq_ID in oligo_dict.keys():
        translated_oligo = str(Seq(oligo_dict[seq_ID]).translate())
        if peptide_dict[seq_ID].replace('X','A').replace('U','C').replace('Z','E').replace('B','D').replace('J','L') != translated_oligo:
            print ('shiiiiiiit')
            print (seq_ID + '\n' + peptide_dict[seq_ID] + '\n' + translated_oligo)


def remove_cutsites(sequence):
    '''Passed a sequence, will return one with  the restriction sites for EcoRI and XhoI replaced/removed with synonymous mutations'''

    #restriction sites
    ecor1 = 'GAATTC'
    xho1 = 'CTCGAG'
    hindIII = 'AAGCCT'
    bser1_f = 'GAGGAG'
    bser1_r = 'CTCCTC'
    mme1_1_f = 'TCCAAC'
    mme1_1_r = 'GTTGGA'
    mme1_2_f = 'TCCGAC'
    mme1_2_r = 'GTCGGA'

    new_sequence = sequence

    for base in range(len(sequence)):
        kmer_6 = sequence[base:base+6]

        if kmer_6 == ecor1:
            #using the base (index) to keep track of the frame as to allow synonymous mutations
            if base%3 == 0: # cut site is in frame, codons read in cut site = GAG TTC, change to GAG TTC
                new_sequence = new_sequence[0:base] + 'GAGTTC'+ new_sequence[base+6:]
            if base%3 == 1: #cut site is shifted one over... meaning codon read is ATT
                new_sequence = new_sequence[0:base] + 'GAATCC'+ new_sequence[base+6:]
            if base%3 == 2:  #cut site is shifted one over... meaning codon read as AAT
                new_sequence = new_sequence[0:base] + 'GAATTC'+ new_sequence[base+6:]

        if kmer_6 == xho1:
            if base%3 == 0: #cut site is in frame, codons read = CTC GAG
                new_sequence = new_sequence[0:base] + 'CTTGAG' + new_sequence[base+6:]
            if base%3 == 1: #cut site is read as XCT CGA GXX. Change CGA to CGT
                new_sequence = new_sequence[0:base] + 'CTCGTG' + new_sequence[base+6:]
            if base%3 == 2: #cut site is read as XXC TCG AGX. Change to AGC (Serine)
                new_sequence = new_sequence[0:base] + 'CAGCAG' + new_sequence[base+6:]

        if kmer_6 == hindIII:
            if base%3 == 0: #cut site is in frame, codons read AAG CCT
                new_sequence = new_sequence[0:base] + 'AAACCA' + new_sequence[base+6:]
            if base%3 == 1: #cut site is read XAA GCC TXX
                new_sequence = new_sequence[0:base] + 'AAGCGT' + new_sequence[base+6:]
            if base%3 == 2: #cut side read as XXA AGC CTX
                new_sequence = new_sequence[0:base] + 'AAGTCT' + new_sequence[base+6:]

    return new_sequence

### variables required for translation ###
bases = ['t', 'c', 'a', 'g']
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))

def translate(seq):
    seq = seq.lower().replace('\n', '').replace(' ', '')
    peptide = ''
    for i in range(0, len(seq), 3):
        codon = seq[i: i+3]
        amino_acid = codon_table.get(codon, '*')
        if amino_acid != '*':
            peptide += amino_acid
        else:
            break
    return peptide



kmer_size = args.k*3
window_size = args.w*3

###Read in the infile###
infile = args.i     #.fasta file
print(infile)

outfile_1 = open(args.o + '.fasta', 'w')
outfile_2 = open(args.o + '.txt', 'w')

###Create dictionary of {seq_ID:peptide_sequence} values
peptide_dict = open_infile(infile)
print ("len(peptide_dict) = " + str(len(peptide_dict.keys())))

### Special Case #1: Just want to split peptide sequences into phage-length seqs, return amino acid seq ###
if(args.p):   #return only a peptide dictionary - do not proceed with script
    cut_peptide_dict = split_kmer_dict(peptide_dict, args.k, args.w)
    for peptide_id in cut_peptide_dict.keys():
        seq_count = 1
        for seq in cut_peptide_dict[peptide_id]:
            outfile_1.write(peptide_id + '_seq' + str(seq_count) + '\n')
            outfile_1.write(seq+ '\n')
            seq_count += 1
    sys.exit()

### Special Case #2: Input is a NT sequence - only do tagging with 5' and 3' sequences and cutsite removal ###
if(args.r):
    nt_dict = peptide_dict   #this is actually NT, so rename for clarity
    for seq_id in nt_dict.keys():
        k = nt_dict[seq_id]
        full_seq = args.fp + k + args.tp
        full_seq_restriction_stripped = remove_cutsites(full_seq)

        if not translate(full_seq) == translate(full_seq_restriction_stripped):
            print("Error: restriction removal caused mismatch at peptide level! Not written to file.")
            print(full_seq)
            print(full_seq_restriction_stripped)

        else: #protein sequence is valid, so write to file
            outfile_1.write(seq_id + '\n')
            outfile_1.write(full_seq_restriction_stripped + '\n')
            outfile_2.write(full_seq_restriction_stripped + '\n\n')
    sys.exit()


### Main / Original case - directly translate AA sequence to NT sequences ###

###Generate a dictionary of {seq_ID:dna_sequence} values
oligo_dict = generate_DNA_dict(peptide_dict)

print ("len(oligo_dict) = " + str(len(oligo_dict.keys())))
sanity_check(oligo_dict,peptide_dict)

###Remove Cutsites
cutsites_removed = {}
for seq_ID in oligo_dict.keys():
    cleaned_seq = remove_cutsites(oligo_dict[seq_ID])
    cutsites_removed[seq_ID] = cleaned_seq

###Make sure that you sanity check!
sanity_check(cutsites_removed, peptide_dict)
print("len(cutsites_removed) = " + str(len(cutsites_removed.keys())))

###Split the full nucleotide sequences into specified kmer sizes
kmers_dict = split_kmer_dict(cutsites_removed, kmer_size, window_size)
print("len(kmers_dict) = " + str(len(kmers_dict.keys())))

###Generate order!
full_seqs = {}
full_seq_restriction_stripped_dict = {}
for seq_ID in kmers_dict.keys():
    kmers = kmers_dict[seq_ID]
    full_seq_restriction_stripped_dict[seq_ID] = []

    kmer_counter = 1  #add kmer id to name (to track ID when original gene has been split into multiple kmers)
    for k in kmers:

        #if not (args.r):
        #    full_seq = args.fp + k + args.tp
        #else:
        #    full_seq = k

        full_seq = args.fp + k + args.tp
        full_seq_restriction_stripped = remove_cutsites(full_seq)
        full_seq_restriction_stripped_dict[seq_ID].append(full_seq_restriction_stripped)
        outfile_1.write(seq_ID + '_seq' + str(kmer_counter) + '\n')
        kmer_counter += 1
        outfile_1.write(full_seq_restriction_stripped + '\n')
        outfile_2.write(full_seq_restriction_stripped + '\n\n')
print ("len(full_seq_restriction_stripped_dict.keys()) = " + str(len(full_seq_restriction_stripped_dict.keys())))

print("completed script.")
