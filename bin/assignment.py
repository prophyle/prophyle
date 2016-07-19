#! /usr/bin/env python3

import os
import shutil
import datetime
import sys
import argparse

from ete3 import Tree

import logging

DEFAULT_FORMAT = 1

def coverage_cigar(hit_list, k):
	qlen=len(hit_list)+k-1

	cov_list=qlen*["X"]
	for i in range(len(hit_list)):
		if hit_list[i]==True:
			for j in range(k):
				cov_list[i+j]="="
	full_cigar="".join(cov_list)
	assert len(full_cigar)==qlen, (qlen, full_cigar)
	x=full_cigar.replace("X=","X\t=").replace("=X","=\tX")
	#print(x,file=sys.stderr)
	parts=x.split()
	y=[]
	l=0
	for z in parts:
		ll=len(z)
		l+=ll
		y.append(str(ll))
		y.append(z[0])
	cigar="".join(y)
	#print(cigar,file=sys.stderr)
	assert l==qlen, (l, qlen, hit_list)
	return cigar

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
		#print("lca=",lca,file=sys.stderr)
		#print(lca.is_root(),lca.children,file=sys.stderr)

		if lca.is_root() and len(lca.children)==1:
			lca=lca.children[0]

		return lca.name

	def name2gi(self,name):
		#print("err",name,file=sys.stderr)
		try:
			gi=self.name_dict[name].gi
		except AttributeError:
			return None
		return gi

	# scores from kmers hits

	"""
		hit1_dict: dict of hit lists - without propagation
		hit2_dict: dict of hit lists - with propagation
	"""
	def assign(self,kmers_assigned_l,simulate_lca=False):
		all_nodes_hit=set()

		hit1_dict=self.dict_from_list(kmers_assigned_l,lca=simulate_lca)
		hit2_dict=hit1_dict.copy()

		for noden in hit1_dict:
			if noden=="0":
				continue

			hit_list=hit1_dict[noden]

			node=self.name_dict[noden]
			while node.up:
				node=node.up
				#print("node up",node.name,file=sys.stderr)
				try:
					hit2_dict[node.name]=hit2_dict[node.name] or hit_list
				except KeyError:
					hit2_dict[node.name]=hit_list

		return hit2_dict

	def print_sam_header(self,file=sys.stdout):
		print("@HD\tVN:1.5\tSO:unsorted",file=file)
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
						rlen=4242,
						as_=as_,
						ur=ur,
						sp=sp,
					),file=file)

	def print_sam_line(self,qname,qlen,rname,krakenmers,hit_list,gi=None,file=sys.stdout):
		flag=0
		pos="1"
		rname2=rname
		#cigar="{}I".format(qlen)
		mapq="60"

		k=qlen-len(hit_list)+1

		if rname is False:
			flag+=4
			rname2="*"
			pos="0"
			cigar="*"
			mapq="0"
		else:
			cigar=coverage_cigar(hit_list,k)

		tags=[]

		score=sum(hit_list)
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

	def print_kraken_line(self,qname,qlen,rname,krakenmers,hit_list,gi=None,file=sys.stdout):
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
		qlen=int(qlen)

		l=[]

		max_hit=None
		gi=None

		blocks=krakenmers.split(" ")
		for b in blocks:
			(ids,count)=b.split(":")
			l.append((ids.split(","),int(count)))

		hit_dict=ti.assign(l,simulate_lca=lca)
		#print(file=sys.stderr)
		#print(a,file=sys.stderr)
		#print(file=sys.stderr)

		#unclassification criterion
		if hit_dict=={}:
			assigned_node=False
		else:
			try:
				del hit_dict["0"]
			except KeyError:
				pass
			max_hit=-1
			noden_m_l=[]
			for noden in hit_dict:
				hit=sum(hit_dict[noden])
				if hit==max_hit:
					noden_m_l.append(noden)
				elif hit>max_hit:
					noden_m_l=[noden]
					max_hit=hit

			#print("final_noden_m_l",noden_m_l,file=sys.stderr)

			if len(noden_m_l)==1:
				#print("non lca",file=sys.stderr)
				assigned_node=noden_m_l[0]
			else:
				#print(noden_m_l,file=sys.stderr)
				assigned_node=ti.lca(noden_m_l)
				#print("lca",assigned_node,file=sys.stderr)
				#print(hit_dict.keys(),file=sys.stderr)
			gi=ti.name2gi(assigned_node)

		#print(hit_dict,file=sys.stderr)
		#print(assigned_node,file=sys.stderr)
		print_line(qname=qname,qlen=qlen,rname=assigned_node,hit_list=hit_dict[assigned_node],krakenmers=krakenmers,gi=gi)
