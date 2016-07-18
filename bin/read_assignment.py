#! /usr/bin/env python3

import os
import shutil
import datetime
import sys
import argparse

from ete3 import Tree

import logging

DEFAULT_FORMAT = 1

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

	def dict_from_list(self,kmers_assigned_l,lca=False):
		d={}

		hit_vector=[]

		npos=sum([x[1] for x in kmers_assigned_l])

		pos=0
		for (noden_l, count) in kmers_assigned_l:			
			if noden_l!=["0"]:
				if lca:
					noden_l=[self.lca(noden_l)]

				v=pos*[False] + count*[True] + (npos-pos-count)*[False]

				for noden in noden_l:
					try:
						assert len(d[noden])==len(v)
						d[noden]=d[noden] or v
					except KeyError:
						d[noden]=v

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
		return lca.name

	def name2gi(self,name):
		#print("err",name,file=sys.stderr)
		try:
			gi=self.name_dict[name].gi
		except AttributeError:
			return None
		return gi

	# scores from kmers hits
	def assign(self,kmers_assigned_l,simulate_lca=False):
		all_nodes_hit=set()

		d=self.dict_from_list(kmers_assigned_l,lca=simulate_lca)
		w=d.copy()

		for noden in d:
			if noden=="0":
				continue
			node=self.name_dict[noden]
			while node.up:
				node=node.up
				if node.name in d:
					w[noden]=w[noden] or d[node.name]
		return w

	def print_sam_header(self,file=sys.stdout):
		print("@HD VN:1.5 SO:unsorted",file=file)
		for node in self.tree.traverse("postorder"):
			self.name_dict[node.name]=node

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
						rlen=1,
						as_=as_,
						ur=ur,
						sp=sp,
					),file=file)

	def print_sam_line(self,qname,qlen,rname,krakenmers,score=None,gi=None,file=sys.stdout):
		flag=0
		pos="1"
		rname2=rname
		cigar="{}I".format(qlen)
		mapq="60"

		if rname is False:
			flag+=4
			rname2="*"
			pos="0"
			cigar="*"
			mapq="0"

		tags=[]
		if score is not None:
			tags.append("AS:i:{}".format(score))

		if gi is not None:
			tags.append("GI:Z:{}".format(gi))

		print("\t".join(
				[
					qname,str(flag),rname2,
					pos,mapq,cigar,
					"*","0", "0","*","*",
					"\t".join(tags)
					#krakenmers
				]
			),file=file)

	def print_kraken_line(self,qname,qlen,rname,krakenmers,score=None,gi=None,file=sys.stdout):
		if rname is False:
			stat="U"
		else:
			stat="C"

		rname2="0" if rname is False else rname

		print("\t".join([stat,qname,rname2,qlen,krakenmers]),file=file)


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

	if form=='sam':
		ti.print_sam_header()

	#print("lets go")


	if form=='kraken':
		print_line=ti.print_kraken_line
	elif form=='sam':
		print_line=ti.print_sam_line


	#ti.process_node(ti.tree.get_tree_root())
	for x in inp_fo:
		x=x.strip()
		stat,qname,_,qlen,krakenmers=x.split("\t")

		l=[]

		max_hit=None
		gi=None

		blocks=krakenmers.split(" ")
		for b in blocks:
			(ids,count)=b.split(":")
			l.append((ids.split(","),int(count)))

		a=ti.assign(l,simulate_lca=lca)
		#print(file=sys.stderr)
		#print(a,file=sys.stderr)
		#print(file=sys.stderr)

		#unclassification criterion
		if a=={}:
			assigned_node=False
		else:
			try:
				del a["0"]
			except KeyError:
				pass
			max_hit=-1
			noden_m_l=[]
			for noden in a:
				hit=sum(a[noden])
				if hit==max_hit:
					noden_m_l.append(noden)
				elif hit>max_hit:
					noden_m_l=[noden]
					max_hit=hit

			if len(noden_m_l)==1:
				assigned_node=noden_m_l[0]
			else:
				#print(noden_m_l,file=sys.stderr)
				assigned_node=ti.lca(noden_m_l)
			gi=ti.name2gi(assigned_node)

		print_line(qname=qname,qlen=qlen,rname=assigned_node,krakenmers=krakenmers,score=max_hit,gi=gi)
