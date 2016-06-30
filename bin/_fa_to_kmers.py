#! /usr/bin/env python3

import sys
import os
import argparse
import re
from itertools import product, repeat
from Bio import SeqIO


reg_splitting=re.compile("[^ACGT]")

comp_dict={
		"A":"T",
		"C":"G",
		"G":"C",
		"T":"A",
	}

def reverse_complement_str(dna):
	reverse_complement="".join([comp_dict[x] for x in dna[::-1]])
	return reverse_complement

def get_kmers_from_fasta(fasta_fn, k, reverse_complement=True):
	kmers=set()

	reg_splitting=re.compile("[^ACGT]")
	set_of_kmers=set()
	fasta_sequences = SeqIO.parse(fasta_fn,'fasta')
	for fasta_seq in fasta_sequences:
		name, sequence = fasta_seq.id, str(fasta_seq.seq).upper()
		sequences_ok=reg_splitting.split(sequence)
		for seq in sequences_ok:
			for i in range(len(seq)-k+1):
				kmer=seq[i:i+k]
				set_of_kmers.add(kmer)
				if reverse_complement:
					kmer_rc=reverse_complement_str(kmer)
					set_of_kmers.add(kmer_rc)

	return set_of_kmers


parser = argparse.ArgumentParser(description='Generate all possible k-mers of given length.')

parser.add_argument(
		'-i','--input',
		help='input fasta file',
		required=True,
	)

parser.add_argument(
		'-k',
		type=int,
		required=True,
		help='k-mer length',
	)

parser.add_argument(
		'-n','--no-revcomp',
		action='store_true',
		help='do not add reverse complements',
	)

parser.add_argument(
		'-f','--format',
		choices=['txt','fq','fa'],
		default='txt',
		help='output format',
	)

args = parser.parse_args()

def print_fq(i,kmer):
	print('@seq'+str(i))
	print(kmer)
	print('+')
	print(args.k*'~')

def print_fa(i,kmer):
	print('>'+str(i))
	print(kmer+'\n')

def print_txt(i,kmer):
	print(kmer)

if args.format=="txt":
	pr=print_txt

elif args.format=="fq":
	pr=print_fq

elif args.format=="fa":
	pr=print_fa

kmers=get_kmers_from_fasta(args.input,args.k,not args.no_revcomp)
kmers=list(kmers)
kmers.sort()

i = 1
for kmer in kmers:
	kmer=''.join(kmer)
	pr(i,kmer)
	i += 1
