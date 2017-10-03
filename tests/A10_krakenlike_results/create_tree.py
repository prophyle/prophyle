from ete3 import Tree

t = Tree('((n1,n2,n3)n4,(n5,n6)n7);', format=1)
t.name='n8'
for node in t:
    node.add_features(fastapath='fastas/'+node.name+'.fa',seqname=node.name+'_seq')
t.write(format=1,features=['fastapath','seqname'],outfile='tree.nw',format_root_node=True)
