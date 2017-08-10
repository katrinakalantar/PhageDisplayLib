
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
#import IPython.html.widgets
#get_ipython().magic(u'matplotlib inline')


# In[81]:

###Functions###
import Bio
from Bio.Seq import Seq


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
        short_list.append(sequence)
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


# In[ ]:




# In[82]:

kmer_size = 62*3                #make this an arg
window_size = 62*3                #make this an arg

###Read in the infile###
infile = 'sequences_greater62.fasta'      #make this an arg

###Create dictionary of {seq_ID:peptide_sequence} values 
peptide_dict = open_infile(infile)
print ("len(peptide_dict) = " + str(len(peptide_dict.keys())))

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

### Add the STREP and FLAG TAGS 
# Strep = Trp-Ser-His-Pro-Gln-Phe-Glu-Lys = WSHPQFEK = tgg agc cat ccg cag ttt gaa aaa  REMOVE FIRST CODON
#FLAG = DYKDDDDK  = gattataaagatgatgatgataag   REMOVE LAST CODON - gat tat aaa gat gat gat gat
strep = 'AGCCATCCGCAGTTCGAGAAA' #Tm = 61.2 salt adjusted
flag = 'GACTACAAGGACGACGATGAT' #Tm =  59.5 salt adjusted 

###Split the full nucleotide sequences into specified kmer sizes
kmers_dict = split_kmer_dict(cutsites_removed, kmer_size, window_size)
print("len(kmers_dict) = " + str(len(kmers_dict.keys())))

###Generate order! 
full_seqs = {}
full_seq_restriction_stripped_dict = {}
outfile_1 = open('agilent_final.fasta', 'w') 
outfile_2 = open('agilent_final.txt', 'w')
for seq_ID in kmers_dict.keys():
    kmers = kmers_dict[seq_ID]
    full_seq_restriction_stripped_dict[seq_ID] = []
    for k in kmers:
        full_seq = strep + k + flag
        full_seq_restriction_stripped = remove_cutsites(full_seq)
        full_seq_restriction_stripped_dict[seq_ID].append(full_seq_restriction_stripped)
        outfile_1.write(seq_ID + '\n')
        outfile_1.write(full_seq_restriction_stripped + '\n')
        outfile_2.write(full_seq_restriction_stripped + '\n\n')
print ("len(full_seq_restriction_stripped_dict.keys()) = " + str(len(full_seq_restriction_stripped_dict.keys())))

print("completed script.")


# In[80]:




# 
