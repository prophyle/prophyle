#! /usr/bin/env python3

import sys
import os
import argparse
import re

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

def load_fasta(fasta_fn):
	print("Loading {}".format(fasta_fn),file=sys.stderr)

	sd={}
	name=None
	seq=[]
	with open(fasta_fn) as f:
		for x in f:
			x=x.upper().strip()
			if x!="":
				if x[0]==">":
					if name!=None:
						sd[name]="".join(seq)
						seq=[]
					name=x[1:]
				else:
					seq.append(x)
	if name!=None:
		sd[name]="".join(seq)

	print("Loaded {} ({} bp)".format(fasta_fn, sum([len(x) for x in sd.values()])),file=sys.stderr)
	return sd

#######
# a = all (both strands)
# f = forward
# r = reversed
# c = canonical
#######
def get_kmers_from_fasta(fasta_fn, k, mode="a"):
	assert mode in ["a","c","r","f"]

	reg_splitting=re.compile("[^ACGT]")
	set_of_kmers=set()

	fasta_dict=load_fasta(fasta_fn)
	print("Extracting k-mers from {} (mode: {})".format(fasta_fn, mode),file=sys.stderr)
	for name in fasta_dict:
		sequence=fasta_dict[name]
		sequences_ok=reg_splitting.split(sequence)
		for seq in sequences_ok:
			if mode=="c":
				for i in range(len(seq)-k+1):
					kmer=seq[i:i+k]
					kmer_rc=reverse_complement_str(kmer)
					set_of_kmers.add(min(kmer,kmer_rc))

			else:
				if mode=="a" or mode=="f":
					for i in range(len(seq)-k+1):
						kmer=seq[i:i+k]
						set_of_kmers.add(kmer)
				if mode=="a" or mode=="r":
					seq_rc=reverse_complement_str(seq)
					for i in range(len(seq_rc)-k+1):
						kmer_rc=seq_rc[i:i+k]
						set_of_kmers.add(kmer_rc)


	print("K-mers extracted from {} (mode: {})".format(fasta_fn, mode),file=sys.stderr)

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
		'-m','--mode',
		choices=["a","c","f","r"],
		default="a",
		help='mode (a = all kmers, c = canonical, f = forward k-mers, r = reversed)',
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

kmers=get_kmers_from_fasta(args.input, args.k, mode=args.mode)
kmers=list(kmers)

print("Sorting all included k-mers from {} ({} kmers)".format(args.input, len(kmers)),file=sys.stderr)
kmers.sort()

i = 1
for kmer in kmers:
	pr(i,kmer)
	i += 1
