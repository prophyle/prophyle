#! /usr/bin/env python3

import os
import shutil
import datetime
import sys
import argparse

from ete3 import Tree

import logging

DEFAULT_FORMAT = 1

class Read:
	def __init__(self, tree, krakline):
		self.tree=tree
		self.load_krakline(krakline)
		self.find_assignments(simulate_lca=False)
		self.annotate_assignments()
		self.print_assignments()

	def load_krakline(self,krakline):
		_,self.qname,_,qlen,self.krakmers=krakline.strip().split("\t")
		self.qlen=int(qlen)
		self.seq=None
		self.qual=None
		self.k=self.qlen+1

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
		self.asgs={}

		for rname in self.hitmasks:
			if rname=="0":
				continue

			self.asgs[rname]={
					'hitmask' : self.hitmasks[rname]
				}

			node=self.tree.name_dict[rname]
			while node.up:
				node=node.up
				#print("node up",node.name,file=sys.stderr)
				try:
					self.asgs[rname]['hitmask']=self.asgs[rname]['hitmask'] or self.hitmasks[node.name]
				except KeyError:
					pass
			self.asgs[rname]['hitmask']="".join(map(str,self.asgs[rname]['hitmask']))



	def annotate_assignments(self):
		"""
		Annotate assignment to a node.

		rname=None => unassigned
		"""

		for rname in self.asgs:
			hitmask=self.asgs[rname]['hitmask']

			"""
			1. hit count
			"""
			self.asgs[rname]['h1']=hitmask.count("1")

			"""
			2. coverage + cigar
			"""
			covlist=self.qlen*[0]
			for i in range(len(hitmask)):
				if hitmask[i]==0:
					for j in range(self.k):
						covlist[i+j]=1
			self.asgs[rname]['covmask']="".join(map(str,covlist))
			self.asgs[rname]['c1']=covlist.count("1")

			x=self.asgs[rname]['covmask'].replace("01","0\t1").replace("10","1\t0")
			y=[]
			for b in x.split():
				y.extend([str(len(b)),"=" if b[0]=="1" else "X"])
			self.asgs[rname]['cigar']="".join(y)

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
		pass

	def print_assignments(self):
		for rname in self.asgs:
			self.print_sam_line(rname)


		#print(file=sys.stderr)
		#print(a,file=sys.stderr)
		#print(file=sys.stderr)

#		#unclassification criterion
#		if hit_dict=={}:
#			assigned_node=False
#			hit_list=None
#			print_line(qname=qname,qlen=qlen,rname=assigned_node,hit_list=hit_list,krakmers=krakmers)
#		else:
#			try:
#				del hit_dict["0"]
#			except KeyError:
#				pass
#			max_hit=-1
#			noden_m_l=[]
#			for noden in hit_dict:
#				hit=sum(hit_dict[noden])
#				if hit==max_hit:
#					noden_m_l.append(noden)
#				elif hit>max_hit:
#					noden_m_l=[noden]
#					max_hit=hit
#
#			#print("final_noden_m_l",noden_m_l,file=sys.stderr)
#
##			if len(noden_m_l)==1:
##				#print("non lca",file=sys.stderr)
##				assigned_node=noden_m_l[0]
##			else:
##				#print(noden_m_l,file=sys.stderr)
##				assigned_node=ti.lca(noden_m_l)
##
##				#print("lca",assigned_node,file=sys.stderr)
##				#print(hit_dict.keys(),file=sys.stderr)
#			for x in noden_m_l:
#				hit_list=hit_dict[x]
#				print_line(qname=qname,qlen=qlen,rname=x,hit_list=hit_list,krakmers=krakmers,)

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

				v=pos*[0] + count*[1] + (npos-pos-count)*[0]

				for noden in noden_l:
					try:
						assert len(d[noden])==len(v)
						d[noden]=d[noden] or v
					except KeyError:
						d[noden]=v

			pos+=count
		return d

#	def load_assignments(self, krakline):
#		x=krakline.strip()
#
#		_,qname,_,qlen,krakenmers=x.split("\t")
#		qlen=int(qlen)
#		l=[]
#
#		max_hit=None
#
#		blocks=krakmers.split(" ")
#		for b in blocks:
#			(ids,count)=b.split(":")
#			l.append((ids.split(","),int(count)))
#
#		hit_dict=ti.assign(l,simulate_lca=lca)
#		#print(file=sys.stderr)
#		#print(a,file=sys.stderr)
#		#print(file=sys.stderr)
#
#		#unclassification criterion
#		if hit_dict=={}:
#			assigned_node=False
#			hit_list=None
#			print_line(qname=qname,qlen=qlen,rname=assigned_node,hit_list=hit_list,krakmers=krakmers)
#		else:
#			try:
#				del hit_dict["0"]
#			except KeyError:
#				pass
#			max_hit=-1
#			noden_m_l=[]
#			for noden in hit_dict:
#				hit=sum(hit_dict[noden])
#				if hit==max_hit:
#					noden_m_l.append(noden)
#				elif hit>max_hit:
#					noden_m_l=[noden]
#					max_hit=hit
#
#			#print("final_noden_m_l",noden_m_l,file=sys.stderr)
#
##			if len(noden_m_l)==1:
##				#print("non lca",file=sys.stderr)
##				assigned_node=noden_m_l[0]
##			else:
##				#print(noden_m_l,file=sys.stderr)
##				assigned_node=ti.lca(noden_m_l)
##
##				#print("lca",assigned_node,file=sys.stderr)
##				#print(hit_dict.keys(),file=sys.stderr)
#			for x in noden_m_l:
#				hit_list=hit_dict[x]
#				print_line(qname=qname,qlen=qlen,rname=x,hit_list=hit_list,krakmers=krakmers,)
#
#		#print(hit_dict,file=sys.stderr)
#	def propagate_assignments(self, asgs):
#		pass
#

#	def krakline2asgs(self, krakline):
#		krakline=krakline.strip()
#
#		b={}
#		_,b['qname'],_,b['qlen'],b['krakmers']=x.split("\t")
#		asg['qlen']=int(asg['qlen'])
#
#		l=[]
#
#		max_hit=None
#
#		blocks=krakmers.split(" ")
#		for b in blocks:
#			(ids,count)=b.split(":")
#			l.append((ids.split(","),int(count)))
#
#		hit_dict=ti.assign(l,simulate_lca=lca)
#		#print(file=sys.stderr)
#		#print(a,file=sys.stderr)
#		#print(file=sys.stderr)
#
#		#unclassification criterion
#		if hit_dict=={}:
#			assigned_node=False
#			hit_list=None
#			print_line(qname=qname,qlen=qlen,rname=assigned_node,hit_list=hit_list,krakmers=krakmers)
#		else:
#			try:
#				del hit_dict["0"]
#			except KeyError:
#				pass
#			max_hit=-1
#			noden_m_l=[]
#			for noden in hit_dict:
#				hit=sum(hit_dict[noden])
#				if hit==max_hit:
#					noden_m_l.append(noden)
#				elif hit>max_hit:
#					noden_m_l=[noden]
#					max_hit=hit
#
#
#	@staticmethod
#	def ann2tags(ann_asg):
#		tags=[]
#
#		try:
#			gi=ann_asg['gi']
#			tags.append("gi:Z:{}".format(gi))
#		except KeyError:
#			pass
#
#		try:
#			taxid=ann_asg['ti']
#			tags.append("ti:Z:{}".format(taxid))
#		except KeyError:
#			pass
#
#		try:
#			h1=ann_asg['h1']
#			tags.append("h1:i:{}".format(h1))
#		except KeyError:
#			pass
#
#		try:
#			h2=ann_asg['h2']
#			tags.append("h2:f:{}".format(h2))
#		except KeyError:
#			pass
#
#		try:
#			c1=ann_asg['c1']
#			tags.append("c1:i:{}".format(c1))
#		except KeyError:
#			pass
#
#		try:
#			c2=ann_asg['c2']
#			tags.append("c2:f:{}".format(c2))
#		except KeyError:
#			pass
#
#		return tags

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


#	"""
#		hit1_dict: dict of hit lists - without propagation
#		hit2_dict: dict of hit lists - with propagation
#	"""
#	def assign(self,kmers_assigned_l,simulate_lca=False):
#		all_nodes_hit=set()
#
#		hit1_dict=self.dict_from_list(kmers_assigned_l,lca=simulate_lca)
#		hit2_dict=hit1_dict.copy()
#
#		for noden in hit1_dict:
#			if noden=="0":
#				continue
#
#			hit2_dict[noden]=hit1_dict[noden]
#
#			node=self.name_dict[noden]
#			while node.up:
#				node=node.up
#				#print("node up",node.name,file=sys.stderr)
#				try:
#					hit2_dict[noden]=hit2_dict[noden] or hit1_list[node.name]
#				except KeyError:
#					pass
#
#		return hit2_dict


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


	def print_kraken_line(self,ann_asg,file=sys.stdout):
		qname=ann_asg['qname']
		qlen=ann_asg['qlen']

		if ann_asg['rname']==False:
			stat="U"
			rname="0"
		else:
			stat="C"
			rname=ann_asg['rname']

		krakmers=ann_asg['krakmers']

		columns=[stat,qname,rname2,qlen,krakmers]

		print("\t".join(columns),file=file)


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

#	if form=='sam':
#		ti.print_sam_header()
#
#
#	if form=='kraken':
#		print_line=ti.print_kraken_line
#	elif form=='sam':
#		print_line=ti.print_sam_line
#
#

	ti.print_sam_header()

	for x in inp_fo:
#		x=x.strip()
#
		read=Read(tree=ti,krakline=x)

#		asg={}
#		_,asg['qname'],_,asg['qlen'],asg['krakmers']=x.split("\t")
#		asg['qlen']=int(asg['qlen'])
#
#		l=[]
#
#		max_hit=None
#
#		blocks=krakmers.split(" ")
#		for b in blocks:
#			(ids,count)=b.split(":")
#			l.append((ids.split(","),int(count)))
#
#		hit_dict=ti.assign(l,simulate_lca=lca)
#		#print(file=sys.stderr)
#		#print(a,file=sys.stderr)
#		#print(file=sys.stderr)
#
#		#unclassification criterion
#		if hit_dict=={}:
#			assigned_node=False
#			hit_list=None
#			print_line(qname=qname,qlen=qlen,rname=assigned_node,hit_list=hit_list,krakmers=krakmers)
#		else:
#			try:
#				del hit_dict["0"]
#			except KeyError:
#				pass
#			max_hit=-1
#			noden_m_l=[]
#			for noden in hit_dict:
#				hit=sum(hit_dict[noden])
#				if hit==max_hit:
#					noden_m_l.append(noden)
#				elif hit>max_hit:
#					noden_m_l=[noden]
#					max_hit=hit
#
#			#print("final_noden_m_l",noden_m_l,file=sys.stderr)
#
##			if len(noden_m_l)==1:
##				#print("non lca",file=sys.stderr)
##				assigned_node=noden_m_l[0]
##			else:
##				#print(noden_m_l,file=sys.stderr)
##				assigned_node=ti.lca(noden_m_l)
##
##				#print("lca",assigned_node,file=sys.stderr)
##				#print(hit_dict.keys(),file=sys.stderr)
#			for x in noden_m_l:
#				hit_list=hit_dict[x]
#				print_line(qname=qname,qlen=qlen,rname=x,hit_list=hit_list,krakmers=krakmers,)
#
#		#print(hit_dict,file=sys.stderr)
#		#print(assigned_node,file=sys.stderr)
#
#		#print_line(qname=qname,qlen=qlen,rname=assigned_node,hit_list=hit_list,krakenmers=krakenmers,gi=gi)
