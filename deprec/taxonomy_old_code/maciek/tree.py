#!/usr/bin/env python
# -*- coding: utf-8 -*-
#!/usr/bin/python -OO -u

##@namespace ithaka.tree
#Taxonomy tree parsing
#@author: "Maciek Sykulski"<macieksk@gmail.com>

import os
import re
#import sys
 
from ithaka.logconf import Config, enum #, data_logger
# create module logger
logger = Config.getLogger('tree')
#logger levels: error,info,warning,debug
#Example log:
#logger.debug("Received from client: {} ClientInfo: {}\nGeoDatagram:{}".format(
#                self.client_address, client_info,
#                self.gdata))#,self.gdata.reserved)

import logging

import ithaka.tree_formatter as trf

import ete3

import csv

from collections import defaultdict
import itertools as it
import multiprocessing
import threading
from ithaka.multiprocmap import parmap

## Class to superpose on a logger to transfer logs by a queue
class LogDummy():
    def __init__(self,q_out):
        self.q_out=q_out
    def __getattr__(self,attr):
        self.attr = attr
        return self        
    def __call__(self,*a,**d):            
        self.q_out.put(('log',logging.__getattribute__(self.attr.upper()),a,d))
        pass            
    pass

def handle_error(logger):
        """Handle an error gracefully.  Overridden to use loger
            The default is to print a traceback and continue.
        """
        try:
            import traceback
            dashln='-'*40
            logger.error('\n{}\nException happened: {}\n{}'.format(
                dashln,traceback.format_exc(),dashln))          
        except Exception as e:
            msg='Exception in handle_error:\n{}'.format(e)
            print(msg)
            logger.error(msg)



def obtain_tree_from_ncbi(db_dir="db",library_dir=Config.library_dir,tree_source="ncbi",
                          nprocs=multiprocessing.cpu_count()):
    try:
        from ete3 import NCBITaxa
        ncbi = NCBITaxa()
        
        os.chdir(db_dir)
        
        id_data_dict = ncbi_taxid_node_data_from_files(library_dir=library_dir,nprocs=nprocs)

        id_tree_f = os.path.join(Config.work_dir,Config.id_tree_f)
        logger.info("Producing Taxonomy subtree from TAXIDS and saving to {}".format(id_tree_f))            
        
        if (len(list(id_data_dict.keys()))>0):
            logger.info("Producing tree topology for {} nodes.".format(len(list(id_data_dict.keys()))))
        else:
            logger.warning("No taxids found! Aborting")
            return        

        if tree_source == "ncbi":   
            tree = ncbi.get_topology(list(id_data_dict.keys()),intermediate_nodes=True)        
        else:
            logger.error("Not implementeted: tree_source={}".format(tree_source))
            return    
        
        extend_tree_with_additional_info(tree,id_data_dict,omit_keys=["taxid"])
        fix_non_leaf_genomic_nodes(tree)
        #Add SUPERROOT since save looses the root name
        parent=type(tree)()
        parent.name = "_SUPERROOT_"
        if (tree.name==""):
            tree.name= "_ROOT_"
        parent.add_child(tree)
        tree = parent
        remove_single_child_nodes(tree)
        save_grepable_newick_format(tree,features=[],outfile=id_tree_f)
    except Exception as e:
        handle_error(logger)
        raise e



def gen_dir_fasta_files(directory,extensions=["fa","ffa","fna","ffn"]):
    logger.info("Begin directory {} walk".format(directory))
    # traverse root directory, and list directories as dirs and files as files
    rcmp = re.compile("({})$".format("|".join(["\."+ext for ext in extensions])))
    for root, dirs, files in os.walk(directory):      
        for file in [f for f in files if re.search(rcmp,f)]:
            yield os.path.join(root,file)		


def ncbi_taxid_node_data_from_files(library_dir = Config.library_dir,
    taxid_map_f=os.path.join(Config.taxonomy_dir,Config.gi_taxid_file),
    nprocs=multiprocessing.cpu_count()):
    
    ##Run thread accepting logs
    #logQueue = multiprocessing.Queue()    
    #def printLogs(logQueue):
        #while True:
            #msg = logQueue.get()
            #if msg is None:
                #break
            #if msg[0]=='log':
                #logger.log(msg[1],*msg[2],**msg[3])
    
    ##Start accepting logs
    #t = threading.Thread(target=printLogs,args=(logQueue,))
    #t.daemon = True
    #t.start()
        
    def parse_fasta(fpath):
        #global logger
        #logger = LogDummy(logQueue)      
        gis = []
        try:        
            logger.info("Parsing fasta file: {}".format(fpath))
            re_description = re.compile("^>")
            with open(fpath,"rt") as f:
                seqnum = 0
                sumlen = 0
                next_addlen = False
                for line in f:
                    #logger.debug("fasta line:{}".format(line))
                    line=line.rstrip()
                    if not re.search(re_description,line): 
                        sumlen+=len(line)
                        continue
                    if next_addlen: gis[-1]['base_len'] = sumlen
                    next_addlen = False
                    sumlen = 0                
                    seqnum += 1
                    #This is info line
                    (seqname,_,_)=line.partition(" ")
                    seqname=seqname[1:] #Omit > character
                    values=seqname.split("|")
                    which_gi = [ i for i,x in enumerate(values) if x == "gi" ]
                    if len(which_gi)<=0: #No gi number
                        logger.warning("Seqid: {} does not have a gi number! Omitting.".format(seqname))
                        continue
                    gi = values[which_gi[0]+1]
                    try: gi = int(gi)
                    except:
                        logger.warning("Seqid: {} does not have an int gi number! Omitting.".format(seqname))
                        continue
                    gis.append({'fastapath':fpath,'infasta_seqnum':seqnum,'seqname':seqname,'gi':gi})
                    next_addlen = True
                if next_addlen: gis[-1]['base_len'] = sumlen
        except:
            handle_error(logger)
        return gis
    
    gis_lists = parmap(parse_fasta,gen_dir_fasta_files(library_dir),nprocs=nprocs)
    #End logging thread
    #logQueue.put(None)
    #t.join()
    
    gis = list(it.chain(*gis_lists))
    logger.info("Done parsing fasta files. Obtained {} records. Sorting by gi numbers...".format(len(gis)))    
    
    gis = sorted(gis,key=lambda x:x['gi'])
    
    #logger.debug("gi.low:{} gi.high:{}".format(gis[0][3],gis[-1][3]))
    
    #Add taxids and save rows tab delimited file
    import csv
    seqid_taxid_map_f = os.path.join(Config.work_dir,Config.seqid_taxid_map_f)
    logger.info("Adding TAXIDS and saving to {}".format(seqid_taxid_map_f))
    id_data_dict=defaultdict(lambda:[])
    with open(seqid_taxid_map_f, 'wt') as csvfile1:
        spamwriter1 = csv.writer(csvfile1, delimiter='\t', lineterminator="\n",
                                quotechar='', quoting=csv.QUOTE_NONE)        
        for i,row in enumerate(map_gi_to_taxid(gis,taxid_map_f)):
            if i==0:
                spamwriter1.writerow(list(row.keys()))
            spamwriter1.writerow(list(row.values()))
            id_data_dict[str(row['taxid'])].append(row)
    
    return id_data_dict
    

def extend_tree_with_additional_info(tree,id_data_dict,omit_keys=["taxid"]):
    cnt = 0
    for n in tree.get_descendants():
        #logger.debug("Node taxid:{}".format(n.taxid))
        if str(n.name) in id_data_dict:
            cnt+=1
            #logger.debug("Found Node taxid:{}".format(n.taxid))
            data_lst = id_data_dict[str(n.name)]
            dct = defaultdict(lambda:[])
            for dl in data_lst:
                for k,v in list(dl.items()):
                    dct[k].append(v)
            for k in list(dct.keys()):
                if k in (["base_len"]+omit_keys): continue
                n.add_feature(k,"@".join(map(str,dct[k])))
            n.add_feature("base_len",sum(dct["base_len"]))
            #n.add_feature("taxid",dct["taxid"][0])
    logger.info("Obtained a tree topology with {} id-mapped nodes.".format(cnt))

def remove_single_child_nodes(tree,preserve_branch_length=True): 
       for n in tree.get_descendants():
            if len(n.children) == 1 and "fastapath" not in n.features:
                n.delete(prevent_nondicotomic=False,
                         preserve_branch_length=preserve_branch_length)

    
def map_gi_to_taxid(gis,taxid_map_f):  
    import gzip
    with gzip.open(taxid_map_f, 'rt') as taxid_map:
        i = 0
        gi = gis[i]['gi']
        for line2 in taxid_map:            
            (gi2, _, taxid2) = line2.partition("\t")
            gi2 = int(gi2)
            #logger.debug("taxid_file|gi:{} taxid:{}".format(gi2,taxid2))             
            while gi <= gi2:
                if gi == gi2:
                    gis[i]['taxid']=int(taxid2.strip())
                    yield gis[i]                    
                else:
                    logger.warning("Could not map gi:{} seqname:'{}' to taxid! Omitting it from the map file. Details:{}".format(gi,gis[i]['seqname'],gis[i]))                    
                i+=1
                if i>=len(gis):
                    return
                if i%500==0:
                    logger.debug("{0:2.1f}% sequence names and gi processed...".format(100.*float(i)/len(gis)))
                gi = gis[i]['gi']                
            

def save_grepable_newick_format(tree, features=None, outfile=None, format=10, is_leaf_fn=None,
              format_root_node=False, dist_formatter=None, support_formatter=None, 
              name_formatter=None):    
    logger.info("Saving tree to filename:{}".format(outfile))
    nw = trf.write_newick(tree, features=features, 
                          format=format,
                          is_leaf_fn=is_leaf_fn,
                          format_root_node=format_root_node, 
                          dist_formatter=dist_formatter,
                          support_formatter=support_formatter, 
                          name_formatter=name_formatter)
    if outfile is not None:
        open(outfile, "wt").write(nw)
    return nw

def fix_non_leaf_genomic_nodes(tree):
    #Obtain the problematic nodes
    inn = [n for n in tree.get_descendants() if len(n.children)>0 and "fastapath" in n.features]
    for node in inn:
        child = node.copy()
        child.children = []        
        node.add_child(child=child)
        node.del_feature("fastapath")
        node.del_feature("seqname")
        node.del_feature("infasta_seqnum")
        node.name = node.name+"_dummyparent"
    #assert(len([n for n in tree.get_descendants() if len(n.children)>0 and "fastapath" in n.features])==0)


def dump_tree(tree_file, outfile,attributes=["est_nkmers"]):
    logger.info("Dumping statistics from taxonomy tree {}".format(tree_file))            
    tree = trf.read_newick(tree_file,format=10)
    with open(outfile, 'wt') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter='\t', lineterminator="\n",
                                quotechar='', quoting=csv.QUOTE_NONE)
        csvwriter.writerow(["name","isleaf"]+attributes)
        for i,node in enumerate(tree.get_descendants()):
            row = [node.name,int(node.is_leaf())]
            row.extend(node.__getattribute__(a) for a in attributes)
            csvwriter.writerow(row)
    


#################### Code in the process of creation


def example_code():
    tree = ncbi.get_topology([9606, 9598, 10090, 7707, 8782],intermediate_nodes=True)
    #print tree.get_ascii(attributes=["sci_name", "rank"])
    tree.write(features=[],outfile="test_tree.newick",format=8)
    save_grepable_newick_format(tree,features=[],outfile="test_tree2.newick")



def mount_memory_fs(directory="."):
    from fs.osfs import OSFS
    from fs.memoryfs import MemoryFS
    from fs.expose import fuse

    home_fs = OSFS(directory)
    home_fs.makedir(Config.ramdrive_dir, allow_recreate=True)
    mp = fuse.mount(MemoryFS(), home_fs.getsyspath(Config.ramdrive_dir))
    return mp

def umount_memory_fs(directory="."):
    from fs.osfs import OSFS
    from fs.expose import fuse

    home_fs = OSFS(directory)    
    fuse.unmount(home_fs.getsyspath(Config.ramdrive_dir))
    