#! /usr/bin/env python3

import os
import shutil
import datetime
import sys
import argparse
import operator


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

	def dict_from_list(self,kmers_assigned_l):
		d={}
		for (noden_l, count) in kmers_assigned_l:
			for noden in noden_l:
				try:
					d[noden]+=count
				except KeyError:
					d[noden]=count
		return d

	def dict_from_list_lca(self,kmers_assigned_l):
		d={}
		for (noden_l, count) in kmers_assigned_l:
			noden=self.lca(noden_l)
			try:
				d[noden]+=count
			except KeyError:
				d[noden]=count
		return d

	def lca(self,noden_l):
		nodes_l=list(map(lambda x:self.name_dict[x],noden_l))
		lca=nodes_l[0].get_common_ancestor(nodes_l)
		return lca.name

	def name2gi(self,name):
		return self.name_dict(name).gi

	def assign(self,kmers_assigned_l,simulate_lca=False):
		all_nodes_hit=set()

		if simulate_lca:
			d=self.dict_from_list_lca(kmers_assigned_l)
		else:
			d=self.dict_from_list(kmers_assigned_l)
		w=d.copy()

		for noden in d:
			node=self.name_dict[noden]
			while node.up:
				node=node.up
				if node.name in d:
					w[noden]+=d[node.name]
		return w

	def print_sam_header(self,file=sys.stdout):
		# @HD VN:1.5 SO:coordinate
		# @SQ SN:ref LN:45
		print("@HD VN:1.5 SO:coordinate",file=file)
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

	def print_sam_line(self,qname,qlen,rname,krakenmers,file=sys.stdout):
		flag=0
		pos="1"
		rname2=rname
		cigar="{}I".format(qlen)
		mapq="1"

		if rname is False:
			flag+=4
			rname2="*"
			pos="0"
			cigar="*"
			mapq="0"

		"r001 99 ref 7 30 8M2I4M1D3M = 37 39 TTAGATAAAGGATACTG *"
		print("\t".join(
				[
					qname,str(flag),rname2,
					pos,mapq,cigar,
					"*","0", "0","*","*"
					#krakenmers
				]
			),file=file)

	def print_kraken_line(self,qname,qlen,rname,krakenmers,file=sys.stdout):
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

	parser.add_argument('-g', '--use-gi',
			action='store_true',
			dest='gi',
			help='output GIs',
		)


	args = parser.parse_args()

	newick_fn=args.newick_fn
	inp_fo=args.input_file
	lca=args.lca
	gi=args.gi
	form=args.format


	ti=TreeIndex(
			tree_newick_fn=newick_fn,
		)

	if form=='sam':
		ti.print_sam_header()

	#print("lets go")


	if form=='kraken':
		printf=ti.print_kraken_line
	elif form=='sam':
		printf=ti.print_sam_line


	#ti.process_node(ti.tree.get_tree_root())
	for x in inp_fo:
		x=x.strip()
		stat,qname,_,qlen,krakenmers=x.split("\t")

		l=[]

		blocks=krakenmers.split(" ")
		for b in blocks:
			(ids,count)=b.split(":")

			if ids=="A" or ids=="0":
				continue

			l.append((ids.split(","),int(count)))

		if l!=[]:

			a=ti.assign(l,simulate_lca=lca)
			stat="C"

			max_hit=-1
			noden_m_l=[]
			for noden in a:
				if a[noden]==max_hit:
					noden_m_l.append(noden)
				elif a[noden]>max_hit:
					noden_m_l=[noden]
					max_hit=a[noden]

			if len(noden_m_l)==1:
				assigned_node=noden_m_l[0]
			else:
				assigned_node=ti.lca(noden_m_l)

			if gi:
				assigned_node=ti.name2gi(assigned_node)
		else:
			assigned_node=False

		printf(qname=qname,qlen=qlen,rname=assigned_node,krakenmers=krakenmers)
		#if form=='kraken':
		#	ti.print_kraken_line(qname=qname,qlen=qlen,rname=assigned_node,krakenmers=krakenmers)
		#elif form=='sam':
		#	ti.print_sam_line(qname=qname,qlen=qlen,rname=assigned_node,krakenmers=krakenmers)
		#print("\t".join([stat,rname,str(assigned_node),rlen,a_kmers]))
