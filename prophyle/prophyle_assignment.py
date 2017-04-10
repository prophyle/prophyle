#! /usr/bin/env python3

"""Prophyle assignment algorithm (a prototype implementation, to be reimplemented in C).

Author:  Karel Brinda <kbrinda@hsph.harvard.edu>

License: MIT

Todo:
	* Allow to use c2 and h2 (when k-mer annotations exist in the NHX tree)
	* Add docstrings
	* Synchronize options with prophyle.py

"""

import os
import shutil
import datetime
import sys
import argparse
import bitarray
import itertools

from ete3 import Tree

DEFAULT_FORMAT = 1

# this should be longer than any possible read
FAKE_CONTIG_LENGTH = 42424242

def cigar_from_mask(mask):
	c=[]
	runs=itertools.groupby(mask)
	for run in runs:
		c.append(str(len(list(run[1]))))
		c.append('=' if run[0] else 'X')
	return "".join(c)


class Read:
	def __init__(self, tree, simulate_lca=False, annotate=False,tie_lca=False,dont_translate_blocks=False):
		self.tree=tree
		self.k=tree.k
		self.simulate_lca=simulate_lca
		self.annotate=annotate
		self.tie_lca=tie_lca
		self.dont_translate_blocks=dont_translate_blocks

	def process_krakline(self,krakline,form,crit):
		self.load_krakline(krakline)
		self.find_assignments()
		#print(self.asgs)
		self.filter_assignments()
		self.print_assignments(form,crit)


	def load_krakline(self,krakline):
		parts=krakline.strip().split("\t")
		self.qname,_,qlen,self.krakmers=parts[1:5]
		self.qlen=int(qlen)
		if len(parts)==7:
			self.seq=parts[5]
			self.qual=parts[6]
		else:
			self.seq="*"
			self.qual="*"

		#self.hitmasks=None
		self.asgs={}

		# list of (count,list of nodes)
		self.kmer_blocks=[]

		b_sum=0
		for b in self.krakmers.split(" "):
			(ids,count)=b.split(":")
			count=int(count)
			b_sum+=count
			rnames=ids.split(",")
			if self.simulate_lca:
				rnames=[self.tree.lca(rnames)]
			self.kmer_blocks.append((rnames,count))
		assert self.qlen==b_sum+self.k-1 or self.qlen<self.k, krakline


	def find_assignments(self):
		# hits before top-down propagation
		hitmasks,covmasks=self.tree.masks_from_kmer_blocks(self.kmer_blocks,simulate_lca=self.simulate_lca)

		rnames=hitmasks.keys()

		# hits after top-down propagation
		for rname in rnames:
			self.asgs[rname]={
					'hitmask' : hitmasks[rname].copy(),
					'covmask' : covmasks[rname].copy(),
				}

			node=self.tree.name_dict[rname]
			for p_rname in self.tree.upnodes_dict[rname] & rnames:
				self.asgs[rname]['hitmask']|=hitmasks[p_rname]
				self.asgs[rname]['covmask']|=covmasks[p_rname]


	def filter_assignments(self):
		"""
		Annotate & filter assignment to a node.

		rname=None => unassigned
		"""

		self.max_hit=0
		self.max_cov=0

		self.max_hit_rnames=[]
		self.max_cov_rnames=[]

		for rname in self.asgs:
			asg=self.asgs[rname]

			"""
			1. hit count
			"""
			hit=asg['hitmask'].count()
			asg['h1']=hit

			if hit>self.max_hit:
				self.max_hit=hit
				self.max_hit_rnames=[rname]
			elif hit==self.max_hit:
				self.max_hit_rnames.append(rname)

			"""
			2. coverage
			"""
			cov=asg['covmask'].count()
			asg['c1']=cov

			if cov>self.max_cov:
				self.max_cov=cov
				self.max_cov_rnames=[rname]
			elif cov==self.max_cov:
				self.max_cov_rnames.append(rname)


	def print_assignments(self, form, crit):
		if crit=="c1":
			winners=self.max_cov_rnames
		elif crit=="h1":
			winners=self.max_hit_rnames

		tie_solved=False

		if self.tie_lca and len(winners)>1:
			tie_solved=True
			lca=self.tree.lca(winners)
			winners=[lca]
			# fix what if this node exists!
			asg=self.asgs[lca]={}
			asg['covcigar']="{}I".format(self.qlen)
			asg['hitcigar']="{}I".format(self.qlen)
			if crit=="h1":
				asg['c1']=-1
				asg['h1']=self.max_hit
			elif crit=="c1":
				asg['c1']=self.max_cov
				asg['h1']=-1

		if len(self.max_hit_rnames)>0:
			for rname in winners:
				asg=self.asgs[rname]
				if form=="sam":
					# compute cigar
					if not tie_solved:
						asg['covcigar']=cigar_from_mask(asg['covmask'])
						asg['hitcigar']=cigar_from_mask(asg['hitmask'])
					self.print_sam_line(rname,self.tree.sam_annotations_dict[rname] if self.annotate else "")
				elif form=="kraken":
					self.print_kraken_line(rname)

		else:
			if form=="sam":
				self.print_sam_line(None)
			elif form=="kraken":
				self.print_kraken_line(None)


	def print_sam_line(self,rname,suffix="",file=sys.stdout):
		tags=[]
		qname=self.qname
		if rname is None:
			flag=4
			rname="*"
			pos="0"
			mapq="0"
			cigar="*"
		else:
			flag=0
			mapq="60"
			pos="1"
			cigar=self.asgs[rname]['covcigar']

		columns=[
				qname,str(flag),rname,
				pos,mapq,cigar,
				"*","0", "0",self.seq,self.qual,
			]

		if rname!="*":
			asg=self.asgs[rname]
			for tag in ['h1','c1']:
				columns.append("".join( [tag,":i:",str(asg[tag])] ))
			columns.append("hc:Z:{}".format(asg['hitcigar']))

		print("\t".join(columns),suffix,file=file,sep="")


	def print_sam_header(self,file=sys.stdout):
		print("@HD\tVN:1.5\tSO:unsorted",file=file)
		for node in self.tree.tree.traverse("postorder"):
			self.tree.name_dict[node.name]=node

			try:
				ur="\tUR:{}".format(node.fastapath)
			except:
				ur=""

			try:
				sp="\tSP:{}".format(node.sci_name)
			except:
				sp=""

			try:
				as_="\tAS:{}".format(node.gi)
			except:
				as_=""

			if node.name!='':
				print("@SQ\tSN:{rname}\tLN:{rlen}{as_}{ur}{sp}".format(
						rname=node.name,
						rlen=FAKE_CONTIG_LENGTH,
						as_=as_,
						ur=ur,
						sp=sp,
					),file=file)


	def print_kraken_line(self,rname,file=sys.stdout):
		if rname is None:
			stat="U"
			rname="0"
		else:
			stat="C"

		if self.simulate_lca:
			lca_rnames=[]
			for [rnames, count] in self.kmer_blocks:
				assert len(rnames)==1

				if rnames[0]=="A" or rnames[0]=="0" or self.dont_translate_blocks:
					taxid=rnames[0]
				else:
					taxid=int(self.tree.taxid_dict[rnames[0]])
				lca_rnames.extend(count*[taxid])
			c=[]
			runs=itertools.groupby(lca_rnames)
			for run in runs:
				c.append("{}:{}".format(
					str(run[0]),
					len(list(run[1]))
				))
			pseudokrakenmers=" ".join(c)

			if stat=="C":
				taxid=str(self.tree.taxid_dict[rname])
			else:
				taxid="0"

			columns=[stat,self.qname,taxid,str(self.qlen),pseudokrakenmers]
			#columns=[stat,self.qname,rname,str(self.qlen)," ".join([ "{}:{}".format(",".join(x[0]),x[1]) for x in self.kmer_blocks])]
		else:
			columns=[stat,self.qname,rname,str(self.qlen),self.krakmers]

		print("\t".join(columns),file=file)


class TreeIndex:

	def __init__(self,tree_newick_fn,k,format=DEFAULT_FORMAT):
		self.tree_newick_fn=tree_newick_fn
		self.tree=Tree(tree_newick_fn,format=1)

		self.k=k

		self.name_dict={}

		self.taxid_dict={}
		self.sam_annotations_dict={}
		self.upnodes_dict={}

		for node in self.tree.traverse("postorder"):
			rname=node.name
			self.name_dict[rname]=node

			# annotations
			tags_parts=[]
			try:
				tags_parts.extend(["\tgi:Z:",node.gi])
			except AttributeError:
				pass

			try:
				tags_parts.extend(["\tti:Z:",node.taxid])
				self.taxid_dict[rname]=node.taxid
			except AttributeError:
				pass

			try:
				tags_parts.extend(["\tsn:Z:",node.sci_name])
			except AttributeError:
				pass

			try:
				tags_parts.extend(["\tra:Z:",node.rank])
			except AttributeError:
				pass

			self.sam_annotations_dict[rname]="".join(tags_parts)

			# set of upper nodes
			self.upnodes_dict[rname]=set()
			while node.up:
				node=node.up
				self.upnodes_dict[rname].add(node.name)


	"""
		kmers_assigned_l:
			list of (list_of_nodes, count)

		dict:
			[
				node => hit_vector,
				node => cov_vector
			]
	"""
	def masks_from_kmer_blocks(self,kmers_assigned_l,simulate_lca=False):
		d_h={}
		d_c={}

		npos=sum([x[1] for x in kmers_assigned_l])

		h_len=npos
		c_len=npos+self.k-1

		pos=0
		for (rname_l, count) in kmers_assigned_l:
			if rname_l!=["0"] and rname_l!=["A"]:
				v_h=bitarray.bitarray(pos*[False] + count*[True] + (npos-pos-count)*[False])
				v_c=bitarray.bitarray(pos*[False] + (count+self.k-1)*[True] + (npos-pos-count)*[False])

				for rname in rname_l:
					try:
						d_h[rname]|=v_h
						d_c[rname]|=v_c
					except KeyError:
						d_h[rname]=v_h.copy()
						d_c[rname]=v_c.copy()

			pos+=count
		return (d_h, d_c)


	def lca(self,noden_l):
		"""Return LCA for a given list of nodes.
		"""
		assert len(noden_l)>0, noden_l
		if len(noden_l)==1:
			return noden_l[0]
		nodes_l=list(map(lambda x:self.name_dict[x],noden_l))
		lca=nodes_l[0].get_common_ancestor(nodes_l)

		if lca.is_root() and len(lca.children)==1:
			lca=lca.children[0]
		assert lca.name!="",[x.name for x in lca.children]
		return lca.name


def assign(
		tree_fn,
		inp_fo,
		lca,
		form,
		k,
		crit,
		annotate,
		tie,
		dont_translate_blocks,
	):
		assert form in ["kraken", "sam"]
		assert k>1
		ti=TreeIndex(
				tree_newick_fn=tree_fn,
				k=k,
			)

		read=Read(
				tree=ti,
				simulate_lca=lca,
				annotate=annotate,
				tie_lca=tie,
				dont_translate_blocks=dont_translate_blocks,
			)
		if form=="sam":
			read.print_sam_header()

		for x in inp_fo:
			read.process_krakline(x,form=form,crit=crit)


def parse_args():
		parser = argparse.ArgumentParser(description='Implementation of assignment algorithm')

		parser.add_argument('-i', '--input',
				type=argparse.FileType('r'),
				required=True,
				dest='input_file',
				help='input file',
			)

		parser.add_argument('-k', '--kmer-size',
				type=int,
				required=True,
				dest='k',
				help='k-mer size',
			)

		parser.add_argument('-n', '--newick-tree',
				type=str,
				metavar='str',
				required=True,
				dest='newick_fn',
				help='newick tree',
			)

		parser.add_argument('-f', '--oformat',
				choices=['kraken','sam'],
				default='kraken',
				dest='format',
				help='format of output',
			)

		parser.add_argument('-m', '--measure',
				choices=['h1','c1'],
				default='h1',
				dest='crit',
				help='measure: h1=hitnumber, c1=coverage',
			)

		parser.add_argument('-l', '--sim-lca',
				action='store_true',
				dest='lca',
				help='simulate LCA',
			)

		parser.add_argument('-a', '--annotate',
				action='store_true',
				dest='annotate',
				help='annotate assignments',
			)

		parser.add_argument('-t', '--tie-lca',
				action='store_true',
				dest='tie',
				help='use LCA when tie (more hits with the same score)',
			)

		parser.add_argument('-d', '--nontransl-blocks',
				action='store_true',
				dest='donttransl',
				help='do not translate blocks from node to tax IDs',
			)
		args = parser.parse_args()
		return args


def main():
		args=parse_args()

		newick_fn=args.newick_fn
		inp_fo=args.input_file
		lca=args.lca
		form=args.format
		k=args.k
		crit=args.crit
		annotate=args.annotate
		tie=args.tie
		d=args.donttransl

		try:
			assign(
					tree_fn=newick_fn,
					inp_fo=inp_fo,
					lca=lca,
					form=form,
					k=k,
					crit=crit,
					annotate=annotate,
					tie=tie,
					dont_translate_blocks=d,
				)
		#Karel: I don't remember why I was considering also IOError here
		#except (BrokenPipeError, IOError):
		except BrokenPipeError:
			# pipe error (e.g., when head is used)
			sys.stderr.close()
			exit(0)


if __name__ == "__main__":
	main()