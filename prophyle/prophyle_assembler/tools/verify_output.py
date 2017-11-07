#! /usr/bin/env python3
"""A program to verify output of prophyle-assembler.

Author: Karel Brinda <kbrinda@hsph.harvard.edu>

Licence: MIT
"""

import sys
import re
from Bio import SeqIO

in1_fn = sys.argv[1]
in2_fn = sys.argv[2]
out1_fn = sys.argv[3]
out2_fn = sys.argv[4]
inter_fn = sys.argv[5]
k = int(sys.argv[6])

reg_splitting = re.compile("[^ACGT]")

comp_dict = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
}


def reverse_complement_str(dna):
    reverse_complement = "".join([comp_dict[x] for x in dna[::-1]])
    return reverse_complement


def get_canonical_kmers_from_fasta(fasta_fn, k):
    kmers = set()

    reg_splitting = re.compile("[^ACGT]")
    set_of_kmers = set()
    fasta_sequences = SeqIO.parse(fasta_fn, 'fasta')
    for fasta_seq in fasta_sequences:
        name, sequence = fasta_seq.id, str(fasta_seq.seq).upper()
        sequences_ok = reg_splitting.split(sequence)
        for seq in sequences_ok:
            for i in range(len(seq) - k + 1):
                kmer = seq[i:i + k]
                kmer_rc = reverse_complement_str(kmer)
                set_of_kmers.add(min(kmer, kmer_rc))
    return set_of_kmers


print("Loading {}".format(in1_fn))
in1 = get_canonical_kmers_from_fasta(in1_fn, k)
print("Loading {}".format(in2_fn))
in2 = get_canonical_kmers_from_fasta(in2_fn, k)
print("Loading {}".format(out1_fn))
out1 = get_canonical_kmers_from_fasta(out1_fn, k)
print("Loading {}".format(out2_fn))
out2 = get_canonical_kmers_from_fasta(out2_fn, k)
print("Loading {}".format(inter_fn))
inter = get_canonical_kmers_from_fasta(inter_fn, k)

print("Is ok", in1 & in2 == inter)

s1 = in1 | in2
s2 = inter | out1 | out2
print("Is ok {} (sizes: {}, {})".format(s1 == s2, len(s1), len(s2)))
print("sym. difference: ", s1 ^ s2)
print()
print("out1 - in1")
print(out1 - in1)
print()
print("out2 - in2")
print(out2 - in2)
