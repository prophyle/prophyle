#! /usr/bin/env python3

import os
import shutil
import datetime
import sys
import argparse
import bitarray
import itertools

from ete3 import Tree



import logging

DEFAULT_FORMAT = 1

class Read:
	def __init__(self, tree, simulate_lca=False, annotate=False):
		self.tree=tree
		self.k=tree.k
		self.simulate_lca=simulate_lca
		self.annotate=annotate

	def process_krakline(self,krakline,form):
		self.load_krakline(krakline)
		self.find_assignments()
		#print(self.asgs)
		self.filter_assignments()
		self.print_assignments(form)


	def load_krakline(self,krakline):
		_,self.qname,_,qlen,self.krakmers=krakline.strip().split("\t")
		self.qlen=int(qlen)
		self.seq=None
		self.qual=None

		#self.hitmasks=None
		self.asgs={}

		# list of (count,list of nodes)
		self.kmer_blocks=[]

		b_sum=0
		for b in self.krakmers.split(" "):
			(ids,count)=b.split(":")
			count=int(count)
			b_sum+=count
			self.kmer_blocks.append((ids.split(","),count))
		assert self.qlen==b_sum+self.k-1, krakline


	def find_assignments(self):
		# hits before top-down propagation
		hitmasks,covmasks=self.tree.masks_from_kmer_blocks(self.kmer_blocks,simulate_lca=self.simulate_lca)
		
		rnames=hitmasks.keys()

		# hits after top-down propagation
		for rname in rnames:
			if rname=="0" or rname=="A":
				continue

			#(hm,cm)=self.hitmasks[rname].copy()

			self.asgs[rname]={
					'hitmask' : hitmasks[rname].copy(),
					'covmask' : covmasks[rname].copy(),
				}

			node=self.tree.name_dict[rname]
			for p_rname in self.tree.upnodes_dict[rname] & rnames:
				self.asgs[rname]['hitmask']|=hitmasks[p_rname]
				self.asgs[rname]['covmask']|=covmasks[p_rname]


			#while node.up:
			#	node=node.up
			#	#print("node up",node.name,file=sys.stderr)
			#	try:
			#		self.asgs[rname]['hitmask']|=hitmasks[node.name]
			#		self.asgs[rname]['covmask']|=covmasks[node.name]
			#	except KeyError:
			#		pass

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
			"""
			1. hit count
			"""
			hit=self.asgs[rname]['hitmask'].count()
			self.asgs[rname]['h1']=hit

			if hit>self.max_hit:
				self.max_hit=hit
				self.max_hit_rnames=[rname]
			elif hit==self.max_hit:
				self.max_hit_rnames.append(rname)

			"""
			2. coverage
			"""
			cov=self.asgs[rname]['covmask'].count()
			self.asgs[rname]['c1']=cov

			if cov>self.max_cov:
				self.max_cov=cov
				self.max_cov_rnames=[rname]
			elif cov==self.max_cov:
				self.max_cov_rnames.append(rname)


			#x=self.asgs[rname]['covmask'].replace("01","0\t1").replace("10","1\t0")
			#y=[]
			#for b in x.split():
			#	y.extend([str(len(b)),"=" if b[0]=="1" else "X"])
			#self.asgs[rname]['cigar']="".join(y)
			#self.asgs[rname]['cigar']="X"


	def annotate_assignment(self, rname, full_annotation=False):
		asg=self.asgs[rname]

		if full_annotation:
			try:
				asg['gi']=self.tree.name_dict[rname].gi
			except AttributeError:
				pass

			try:
				asg['ti']=self.tree.name_dict[rname].taxid
			except AttributeError:
				pass

			try:
				asg['sn']=self.tree.name_dict[rname].sci_name
			except AttributeError:
				pass

			try:
				asg['ra']=self.tree.name_dict[rname].rank
			except AttributeError:
				pass

		c=[]
		runs=itertools.groupby(asg['covmask'])
		for run in runs:
			c.append(str(len(list(run[1]))))
			c.append('=' if run[0] else 'X')
		self.asgs[rname]['cigar']="".join(c)


	def print_assignments(self, form):
		if len(self.max_hit_rnames)>0:
			for rname in self.max_hit_rnames:
				if form=="sam":
					self.annotate_assignment(rname,full_annotation=self.annotate)
					self.print_sam_line(rname)
				elif form=="kraken":
					self.print_kraken_line(rname)

		else:
			if form=="sam":
				self.print_sam_line(None)
			elif form=="kraken":
				self.print_kraken_line(None)


	def print_sam_line(self,rname,file=sys.stdout):
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
			cigar=self.asgs[rname]['cigar']

		columns=[
				qname,str(flag),rname,
				pos,mapq,cigar,
				"*","0", "0","*","*",
			]
		columns.extend(self.sam_tags(rname))
		print("\t".join(columns),file=file)


	def sam_tags(self,rname):
		tags=[]

		if rname!="*":
			asg=self.asgs[rname]
			if self.annotate:
				for tag in ['gi','ti','sn','ra']:
					try:
						tags.append("".join( [tag,":Z:",asg[tag]] ))
					except KeyError:
						pass


			for tag in ['h1','c1']:
				tags.append("".join( [tag,":i:",str(asg[tag])] ))

		return tags

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
						rlen=4242,
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
		columns=[stat,self.qname,rname,str(self.qlen),self.krakmers]
		print("\t".join(columns),file=file)

class TreeIndex:

	def __init__(self,tree_newick_fn,k,format=DEFAULT_FORMAT):
		self.tree_newick_fn=tree_newick_fn
		self.tree=Tree(tree_newick_fn,format=1)

		self.k=k

		self.name_dict={}
		self.upnodes_dict={}

		for node in self.tree.traverse("postorder"):
			rname=node.name
			self.name_dict[rname]=node
			self.upnodes_dict[rname]=set()
			while node.up:
				node=node.up
				self.upnodes_dict[rname].add(node.name)


		#print (self.name_dict)

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
		c_len=npos+k-1

		pos=0
		for (rname_l, count) in kmers_assigned_l:
			if rname_l!=["0"] and rname_l!=["A"]:
				if lca:
					rname_l=[self.lca(rname_l)]

				v_h=bitarray.bitarray(pos*[False] + count*[True] + (npos-pos-count)*[False])
				v_c=bitarray.bitarray(pos*[False] + (count+k-1)*[True] + (npos-pos-count)*[False])

				#assert len(v_h)==h_len, v_h
				#assert len(v_c)==c_len, v_c

				for rname in rname_l:
					try:
						#assert len(d_h[noden])==h_len
						#assert len(d_c[noden])==c_len
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

		return lca.name



if __name__ == "__main__":

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

	args = parser.parse_args()

	newick_fn=args.newick_fn
	inp_fo=args.input_file
	lca=args.lca
	form=args.format
	k=args.k
	annotate=args.annotate


	ti=TreeIndex(
			tree_newick_fn=newick_fn,
			k=k,
		)

	read=Read(
			tree=ti,
			simulate_lca=lca,
			annotate=annotate
		)
	if form=="sam":
		read.print_sam_header()
	for x in inp_fo:
		read.process_krakline(x,form=form)

