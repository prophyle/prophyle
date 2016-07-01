#! /usr/bin/env python3

import sys, os, argparse, re
from collections import deque
from ete3 import PhyloTree, NCBITaxa

def read_buf(file, buf_size=536870912):
	buf = file.read(buf_size)+file.readline()
	while buf:
		yield buf
		buf = file.read(buf_size)+file.readline()

def index_of(taxid, taxa_list):
	for i, seq in enumerate(taxa_list):
		if taxid == seq[5]:
			return i
	return -1

parser = argparse.ArgumentParser(
	description='Assigns taxids to the sequences in the fasta indexes contained '
				'in the input directory and its subdirectories, and builds a '
				'taxonomic tree in the New Hampshire 0 newick format')

parser.add_argument('library_dir', help = 'directory containing the .fai files')
parser.add_argument('taxid_map_f', help = 'map of gis to taxid from NCBI database')
parser.add_argument('-o', '--output-file',
					type=str,
					metavar='output_file',
					default = 'taxonomic_tree.nw',
					dest='output_file',
					help = 'output file (default: taxonomic_tree.nw)')
parser.add_argument('-e', '--error-file',
					type=str,
					metavar='error_file',
					default = 'error_tree.log',
					dest='error_file',
					help = 'error log file (default: error_tree.log)')
parser.add_argument('-a', '--assign-only',
					type=str,
					metavar='dmp_file',
					default = None,
					dest='dmp_file',
					help = 'only assign the taxids to the sequences and store them in the given file')
parser.add_argument('-b', '--build_tree',
					type=str,
					metavar='assignments_file',
					default = None,
					dest='assignments_file',
					help = 'build the tree from the given file (previous output of --assign-only)')

args = parser.parse_args()
library_dir = args.library_dir
if not library_dir.endswith('/'):
	library_dir += '/'
taxid_map_f = args.taxid_map_f
output_file = args.output_file
error_file = args.error_file
dmp_file = args.dmp_file
assignments_file = args.assignments_file

error = open(error_file, 'w')
ass_seqs = []

if assignments_file is not None:
	with open(assignments_file, 'r') as seq_taxid:
	    for line in seq_taxid:
	        values = line.split("\t")
	        ass_seqs.append([str(values[0]),str(values[1]),str(values[2]),str(values[3]),str(values[4]),int(values[5])])
else:
	seqs = []
	skipped = 0
	for dirpath, dirnames, filenames in os.walk(library_dir):
		for filename in (f for f in filenames if f.endswith(".fai")):
			fn = os.path.join(dirpath, filename)
			with open(fn, 'r') as faidx:
				for seq in faidx:
					values = seq.split("\t")
					seqname, seqlenght, offset, _, _ = values
					split_seqname = seqname.split("|")
					gi = split_seqname[2]
					try:
						seqs.append([fn[fn.find("library"):],
									str(seqname), str(seqlenght),
									str(offset), int(gi)])
					except:
						if skipped == 0:
							error.write("NOT ACQUIRED:\n")
						error.write(fn + " " + str(seqname) + "\n")
						skipped += 1
						pass

	seqs = sorted(seqs, key = lambda x:x[4], reverse = True)
	seqs_no = len(seqs)
	print("Acquired " + str(seqs_no) + " seqs (" + str(skipped) + " skipped)")

	skipped = 0
	with open(taxid_map_f, 'r') as taxid_map:
		found = True
		finished = False
		seq = [-1,-1,-1,-1,-1]
		for buf in read_buf(taxid_map):
			if finished: break
			for line in buf.splitlines():
				(gi, _, taxid) = line.partition("\t")
				gi = int(gi)
				while seqs and seq[4] < gi:
					if not found:
						if skipped == 0:
							error.write("NOT ASSIGNED:\n")
						error.write(str(seq[0]) + " " + str(seq[1])+"\n")
						skipped += 1
					else:
						found = False
					seq = seqs.pop()
				if seq[4] == gi:
					found = True
					ass_seqs.append(seq.append(taxid.strip()))
				if not seqs and seq[4] < gi:
					if not found:
						if skipped == 0:
							error.write("\n\nNOT ASSIGNED:\n\n")
						error.write(str(seq[0]) + " " + str(seq[1])+"\n")
						skipped += 1
					finished = True
					break

	print("Assigned " + str(seqs_no-skipped) + " seqs (" + str(skipped) + " skipped)")

	if dmp_file is not None:
		with open(dmp_file, 'w') as output:
			for seq in ass_seqs:
				print(seq)
				output.write("\t".join(map(str,seq))+"\n")
		print("Assignments written to " + dmp_file)
		print("Launch again with --build_tree " + dmp_file +
				" to build a taxonomic tree from them")
		error.close()
		sys.exit(0)

ass_seqs = sorted(ass_seqs, key = lambda x:x[5])
taxids = [-1]
prec = 0
for seq in ass_seqs:
	tid = seq[5]
	if tid != taxids[prec]:
		taxids.append(tid)
		prec += 1
del taxids[0]

ncbi = NCBITaxa()
topo = ncbi.get_topology(taxids)
new_id = 1
count = 0
library_dir = library_dir[library_dir.find("library"):]
for node in topo.traverse("preorder"):
	node.name = new_id
	new_id += 1
	i = index_of(node.taxid, ass_seqs)
	if i != -1:
		gi = ""
		seqname = ""
		fastapath = ""
		infasta_seqnum = ""
		while i < len(ass_seqs) and ass_seqs[i][5] == node.taxid:
			if gi == "":
				fastapath = ass_seqs[i][0]
				seqname = ass_seqs[i][1]
				base_len = ass_seqs[i][2]
				infasta_offset = ass_seqs[i][3]
				gi = str(ass_seqs[i][4])
			else:
				fastapath += "@"+ass_seqs[i][0]
				seqname += "@"+ass_seqs[i][1]
				base_len += "@"+ass_seqs[i][2]
				infasta_offset += "@"+ass_seqs[i][3]
				gi += "@"+str(ass_seqs[i][4])
			i += 1
			count += 1
		node.add_features(fastapath = fastapath, seqname = seqname, base_len = base_len,
							infasta_offset = infasta_offset, gi = gi)

print("Built taxonomic tree for " + str(count) + " sequences")

topo.write(features=["lineage", "named_lineage", "seqname", "dist", "name",
					"support", "taxid", "rank", "base_len", "fastapath",
					"sci_name", "common_name", "infasta_offset", "gi"],
			outfile=output_file)

error.close()
