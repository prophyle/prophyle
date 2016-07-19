#! /usr/bin/env python3

import os
import shutil
import datetime
import sys
import argparse
import numpy as np

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
		self.annotate_assignments()
		self.filter_assignments()
		self.print_assignments()


	def load_krakline(self,krakline):
		_,self.qname,_,qlen,self.krakmers=krakline.strip().split("\t")
		self.qlen=int(qlen)
		self.seq=None
		self.qual=None
		self.k=self.qlen+1

		self.hitmasks=None
		self.asgs={}

		# list of (count,list of nodes)
		self.kmer_blocks=[]
		for b in self.krakmers.split(" "):
			(ids,count)=b.split(":")
			count=int(count)
			self.k-=count
			self.kmer_blocks.append((ids.split(","),count))


	def find_assignments(self,simulate_lca):
		# hits before top-down propagation
		self.hitmasks=self.tree.hitmasks_from_kmer_blocks(self.kmer_blocks,lca=simulate_lca)
		# hits after top-down propagation

		for rname in self.hitmasks:
			if rname=="0":
				continue

			self.asgs[rname]={
					'hitmask' : self.hitmasks[rname].copy()
				}

			node=self.tree.name_dict[rname]
			while node.up:
				node=node.up
				#print("node up",node.name,file=sys.stderr)
				try:
					self.asgs[rname]['hitmask']+=self.hitmasks[node.name]
				except KeyError:
					pass

	def annotate_assignments(self):
		"""
		Annotate assignment to a node.

		rname=None => unassigned
		"""

		for rname in self.asgs:
			"""
			1. hit count
			"""
			self.asgs[rname]['h1']=sum(list(self.asgs[rname]['hitmask']))

			"""
			2. coverage + cigar
			"""
			kmask=np.array(self.k*[True])
			self.asgs[rname]['covmask']=np.convolve(self.asgs[rname]['hitmask'],kmask)
			#print(self.asgs[rname]['hitmask'])
			#print(kmask)
			#print(self.asgs[rname]['covmask'])
			self.asgs[rname]['c1']=sum(list(self.asgs[rname]['covmask']))

			#x=self.asgs[rname]['covmask'].replace("01","0\t1").replace("10","1\t0")
			#y=[]
			#for b in x.split():
			#	y.extend([str(len(b)),"=" if b[0]=="1" else "X"])
			#self.asgs[rname]['cigar']="".join(y)
			self.asgs[rname]['cigar']="X"

			"""
			3. assign taxid & gi
			"""
			try:
				self.asgs[rname]['gi']=self.name_dict[rname].gi
			except AttributeError:
				pass

			try:
				self.asgs[rname]['ti']=self.name_dict[rname].taxid
			except AttributeError:
				pass


	def filter_assignments(self):
		hit_counts=[self.asgs[x]['h1'] for x in self.asgs]
		if len(hit_counts)>0:
			self.max_hit=max(hit_counts)
		else:
			self.max_hit=None


	def print_assignments(self):
		if self.max_hit is not None:			
			for rname in self.asgs:
				if self.asgs[rname]['h1']==self.max_hit:
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
			node => hit_vector
	"""
	def hitmasks_from_kmer_blocks(self,kmers_assigned_l,lca=False):
		d={}

		hit_vector=[]

		npos=sum([x[1] for x in kmers_assigned_l])

		pos=0
		for (noden_l, count) in kmers_assigned_l:			
			if noden_l!=["0"]:
				if lca:
					noden_l=[self.lca(noden_l)]

				v=np.array(pos*[False] + count*[True] + (npos-pos-count)*[False])

				for noden in noden_l:
					try:
						assert len(d[noden])==len(v)
						d[noden]+=v
					except KeyError:
						d[noden]=v.copy()

			pos+=count
		return d


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

