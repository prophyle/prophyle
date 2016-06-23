#! /usr/bin/env python3

class Node:
	def __init__(self,id,parent,children=(),data={}):
		self.id=id
		self.parent=parent
		self.children=children
		self.data=data

	def repr(self):
		return "{}:({})".format(self.id,",".join([x.repr() for x in self.children]))

class TaxonomyTree:
	def __init__(self):
		# node: [id, taxid, taxid_parent, taxid_children]
		self.root = None
		self.nodes = []
		self.nodes_dict = {}
		self.leaves = []

		#self.ranks = {}
		#self.seqnames = {}
		#self.gis = {}

	def get_root(self):
		return self.root

	def save_tree(self,filename):
		pass

	def load_tree(self,filename):
		pass
	
	def repr(self):
		if self.root is not None:
			return self.root.repr()

	# from file nodes.tmp
	def load_from_ncbi_dmp(self,taxid_dmp):
		self.root = None
		self.nodes_dict = {}
		self.leaves = []

		# {tad_id : ( data, children )}
		records = {}
		with open(taxid_dmp) as dmp:
			for line in dmp:
				values=line.split("|")
				tax_id=int(values[0])
				parent_tax_id=int(values[1])
				rank=values[2].strip()
				# parent's record
				try:
					# child's id
					records[parent_tax_id]=\
						(
							records[parent_tax_id][0],
							records[parent_tax_id][1]+(tax_id,),
						)
				except KeyError:
					# new parent's record
					records[parent_tax_id]=\
						(
							(),
							(tax_id,),
						)
				# own record
				try:
					# only data
					records[tax_id]=\
						(
							(tax_id,parent_tax_id,rank),
							records[tax_id][1],
						)
				except KeyError:
					# complete record
					records[tax_id]=\
						(
							(tax_id,parent_tax_id,rank),
							(),
						)

		print("dmp loaded")

		# correction for root
		root_i=records[1][0]
		root_children=records[1][1]
		records[1]=\
			(
				(root_i[0],None,""),
				[x for x in root_children if x!=1],
			)

		#print(records)

		def _add_tax_id_subbranch(tax_id,parent_id):
			assert (parent_id is None) == (self.root is None)
			#print("adding subbranch {} {}".format(tax_id,parent_id))
			#print(records[tax_id])
			((tax_id,parent_tax_id,rank),children_ids) = records[tax_id]
			node = Node(
						id=tax_id,
						parent=self.get_node_from_id(parent_id),
					)
			self.add_node(node)
			#print("children ids",children_ids)
			for child_id in children_ids:
				_add_tax_id_subbranch(child_id,tax_id)

		_add_tax_id_subbranch(1,None)


	def add_node(self, node):
		assert node.id not in self.nodes_dict.keys()
		assert len(node.children)==0

		if self.root is None:
			assert node.parent is None
			self.root=node
		else:
			assert node.parent in self.nodes, "each new node must have a parent"
			try:
				self.leaves.remove(parent)
			except:
				pass
			node.parent.children = node.parent.children + (node,)

		self.nodes.append(node)
		self.nodes_dict[node.id]=node
		self.leaves.append(node)

	def get_node_from_id(self,id):
		if id is None:
			return None
		else:
			return self.nodes_dict[id]

	def remove_node(self,node):
		for child in node.children:
			self.remove_node(child)
		self.leaves.append(node.parent)

		del self.nodes_dict[node.id]
		self.nodes.remove(node)
		self.leaves.remove(node)

	def split_node(self,node,children1,children2):
		pass

	def remove_leaf_path(self,leaf_id):
		assert(self.is_leaf(id))
		self.remove_node(leaf_id)
		while(leaf_id!=0):
			pass

	def is_leaf(self,node_id):
		return node_id in self.leaves

	def update_node(self,seqname=None,gid=None):
		pass



tt = TaxonomyTree()
tt.load_from_ncbi_dmp("nodes.dmp")
print("tree loaded")
print(tt.repr())
