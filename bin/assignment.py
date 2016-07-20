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
	def __init__(self, tree):
		self.tree=tree

	def process_krakline(self,krakline):
		self.load_krakline(krakline)
		self.find_assignments(simulate_lca=False)
		#print(self.asgs)
		self.filter_assignments()
		self.print_assignments()


	def load_krakline(self,krakline):
		_,self.qname,_,qlen,self.krakmers=krakline.strip().split("\t")
		self.qlen=int(qlen)
		self.seq=None
		self.qual=None
		self.k=self.qlen+1

		#self.hitmasks=None
		self.asgs={}

		# list of (count,list of nodes)
		self.kmer_blocks=[]

		b_sum=0
		for b in self.krakmers.split(" "):
			(ids,count)=b.split(":")
			count=int(count)
			self.k-=count
			b_sum+=count
			self.kmer_blocks.append((ids.split(","),count))
		assert self.qlen==b_sum, krakline


	def find_assignments(self,simulate_lca):
		# hits before top-down propagation
		hitmasks,covmasks=self.tree.masks_from_kmer_blocks(self.kmer_blocks,lca=simulate_lca,k=self.k)
		
		# hits after top-down propagation
		for rname in hitmasks:
			if rname=="0" or rname=="A":
				continue

			#(hm,cm)=self.hitmasks[rname].copy()

			self.asgs[rname]={
					'hitmask' : hitmasks[rname].copy(),
					'covmask' : covmasks[rname].copy(),
				}

			node=self.tree.name_dict[rname]
			while node.up:
				node=node.up
				#print("node up",node.name,file=sys.stderr)
				try:
					self.asgs[rname]['hitmask']|=hitmasks[node.name]
					self.asgs[rname]['covmask']|=covmasks[node.name]
				except KeyError:
					pass

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


	def annotate_assignment(self, rname):
		try:
			self.asgs[rname]['gi']=self.tree.name_dict[rname].gi
		except AttributeError:
			pass

		try:
			self.asgs[rname]['ti']=self.tree.name_dict[rname].taxid
		except AttributeError:
			pass

		c=[]
		runs=itertools.groupby(self.asgs[rname]['covmask'])
		for run in runs:
			c.append(str(len(list(run[1]))))
			c.append('=' if run[0] else 'X')
		self.asgs[rname]['cigar']="".join(c)


	def print_assignments(self):
		if len(self.max_hit_rnames)>0:
			for rname in self.max_hit_rnames:
				self.annotate_assignment(rname)
				self.print_sam_line(rname)
		else:
			self.print_sam_line(None)


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

		try:
			gi=self.asgs[rname]['gi']
			tags.append("gi:Z:{}".format(gi))
		except KeyError:
			pass

		try:
			taxid=self.asgs[rname]['ti']
			tags.append("ti:Z:{}".format(taxid))
		except KeyError:
			pass

		try:
			h1=self.asgs[rname]['h1']
			tags.append("h1:i:{}".format(h1))
		except KeyError:
			pass

		try:
			h2=self.asgs[rname]['h2']
			tags.append("h2:f:{}".format(h2))
		except KeyError:
			pass

		try:
			c1=self.asgs[rname]['c1']
			tags.append("c1:i:{}".format(c1))
		except KeyError:
			pass

		try:
			c2=self.asgs[rname]['c2']
			tags.append("c2:f:{}".format(c2))
		except KeyError:
			pass

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


	def print_kraken_line(self,ann_asg,file=sys.stdout):
		if rname is None:
			stat="U"
			rname="0"
		else:
			stat="C"
		columns=[stat,self.qname,rname,self.qlen,self.krakmers]
		print("\t".join(columns),file=file)

class TreeIndex:

	def __init__(self,tree_newick_fn,format=DEFAULT_FORMAT):
		self.tree_newick_fn=tree_newick_fn
		self.tree=Tree(tree_newick_fn,format=1)

		self.name_dict={}

		for node in self.tree.traverse("postorder"):
			self.name_dict[node.name]=node

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
	def masks_from_kmer_blocks(self,kmers_assigned_l,k,lca=False):
		d_h={}
		d_c={}

		assert k>1

		npos=sum([x[1] for x in kmers_assigned_l])

		h_len=npos
		c_len=npos+k-1

		pos=0
		for (noden_l, count) in kmers_assigned_l:
			print("block",noden_l,count)
			if noden_l!=["0"] and noden_l!=["A"]:
				if lca:
					noden_l=[self.lca(noden_l)]

				v_h=bitarray.bitarray(pos*[False] + count*[True] + (npos-pos-count)*[False])
				v_c=bitarray.bitarray(pos*[False] + (count+k-1)*[True] + (npos-pos-count)*[False])

				assert len(v_h)==h_len
				print(k,(count+k-1))
				assert len(v_c)==c_len, v_c

				for noden in noden_l:
					try:
						assert len(d_h[noden])==h_len
						assert len(d_c[noden])==c_len
						d_h[noden]|=v_h
						d_c[noden]|=v_c
					except KeyError:
						d_h[noden]=v_h.copy()
						d_c[noden]=v_c.copy()

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

	args = parser.parse_args()

	newick_fn=args.newick_fn
	inp_fo=args.input_file
	lca=args.lca
	form=args.format


	ti=TreeIndex(
			tree_newick_fn=newick_fn,
		)

	read=Read(tree=ti)
	read.print_sam_header()
	for x in inp_fo:
		read.process_krakline(x)

