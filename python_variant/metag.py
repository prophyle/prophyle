#! /usr/bin/env python3

import os
import re
from Bio import SeqIO

comp_dict={
		"A":"T",
		"C":"G",
		"G":"C",
		"T":"A",	
	}
	
dna_alphabet=set(["A","C","G","T"])

def reverse_complement(dna):
	reverse=dna[::-1]
	reverse_complement="".join([comp_dict[x] for x in reverse])
	return reverse_complement
	
def canonical_repr(dna):
	assert isinstance(dna,str)
	forward=dna[::-1]
	reverse_complement=reverse_complement(forward)
	return forward if forward < reverse_complement else reverse_complement

reg_dna=re.compile(r"^[ACGT]+$")
reg_seed=re.compile(r"^[#-]+$")
def is_dna(dna):
	#return set([x for x in dna]).issubset(dna_alphabet)
	return bool(reg_dna.match(dna))

def set_from_fasta(fasta_fn,k,canonical=False):
	reg_splitting=re.compile("[^ACGT]")
	set_of_kmers=set()
	fasta_sequences = SeqIO.parse(open(fasta_fn),'fasta')
	for fasta_seq in fasta_sequences:
		name, sequence = fasta_seq.id, str(fasta_seq.seq)
		sequences_ok=reg_splitting.split(sequence)
		for seq in sequences_ok:
			for i in range(len(seq)-k+1):
				kmer=seq[i:i+k]
				if canonical:
					set_of_kmers.add(canonical_repr(kmer))
				else:
					set_of_kmers.add(kmer)
					set_of_kmers.add(reverse_complement(kmer))
	return set_of_kmers


def set_from_fasta_deprec(fasta_fn,k,canonical=False):
	set_of_kmers=set()
	fasta_sequences = SeqIO.parse(open(fasta_fn),'fasta')
	for fasta in fasta_sequences:
		name, sequence = fasta.id, str(fasta.seq)
		for i in range(len(sequence)-k+1):
			kmer=sequence[i:i+k]
			if is_dna(kmer):
				if canonical:
					set_of_kmers.add(canonical_repr(kmer))
				else:
					set_of_kmers.add(kmer)
					set_of_kmers.add(reverse_complement(kmer))
	return set_of_kmers


def set_from_fasta_n(fasta_fn,k):
	dna_dict={
			"A": 0,
			"C": 1,
			"G": 2,
			"T": 3,
		}

	import numpy

	with open(fasta_fn) as f:
		lines = [line.strip().upper() for line in f if line[0]!=">"]
		fasta = "".join(lines)
		lines=[]
		nseed=numpy.array([4**i for i in range(k)])
		index=set()
		for i in range(len(fasta)-len(nseed)+1):
			substr=fasta[i:i+len(nseed)]
			# valid dna string?
			if is_dna(substr):
				nkmer = numpy.array( [ dna_dict[x] for x in substr ] )
				number = numpy.dot( nseed, nkmer )
				index.add(number)
	return index



def set_to_fasta(fasta_fn,set_of_kmers,assemble=False):
	with open(fasta_fn,"w+") as fasta_fo:
		if assemble:
			set_of_kmers|=set([reverse_complement(kmer) for kmer in set_of_kmers])
			#print(set_of_kmers)
			i=0
			while len(set_of_kmers)>0:
				initial_kmer=min(set_of_kmers)
				initial_kmer_rc=reverse_complement(initial_kmer)
				k=len(initial_kmer)
				set_of_kmers.remove(initial_kmer)
				if initial_kmer!=initial_kmer_rc:
					set_of_kmers.remove(initial_kmer_rc)
				contig=list(initial_kmer)
				extending=True
				while extending:
					extending=False
					for right_ext in ["A","C","G","T"]:
						cand_kmer="".join(contig[-k+1:]+[right_ext])
						if cand_kmer in set_of_kmers:
							contig.append(right_ext)
							set_of_kmers.remove(cand_kmer)					
							cand_kmer_rc=reverse_complement(cand_kmer)
							if cand_kmer!=cand_kmer_rc:
								set_of_kmers.remove(cand_kmer_rc)
							extending=True
							break
				fasta_fo.write("".join([">contig_",fasta_fn,"_",str(i),os.linesep,"".join(contig),os.linesep]))
				i+=1				
		else:
			for kmer,i in enumerate(set_of_kmers):
				fasta_fo.write("".join([">contig",str(i),os.linesep,kmer,os.linesep]))


def and_diff(source_fasta_1_fn, source_fasta_2_fn, filtered_fasta_1_fn, filtered_fasta_2_fn, upper_fasta_fn):
	set_1=set_from_fasta(source_fasta_1_fn)
	set_2=set_from_fasta(source_fasta_2_fn)
	upper_set=set_1&set_2
	set_1-=upper_set
	set_2-=upper_set
	set_to_fasta(filtered_fasta_1_fn,set_1)
	set_to_fasta(filtered_fasta_2_fn,set_2)
	set_to_fasta(upper_fasta_fn,upper_set)

if __name__ == "__main__":
	fasta_fn="../tests/Borrelia_garinii.fa"
	k=2
	kmers=set_from_fasta(fasta_fn,k=2)
	set_to_fasta("test.fa",kmers,assemble=True)