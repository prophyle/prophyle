# ProPhyle - a metagenomic classifier

[![Build Status](https://travis-ci.com/karel-brinda/ProPhyle.svg?token=LzzDiQkWWqF4hBjZahmQ&branch=master)](https://travis-ci.com/karel-brinda/ProPhyle)

## Getting started

### Prerequisities

* GCC 4.8+
* ZLib
* Python 3 with ete3 library
* SamTools

#### Recommended way of installation using [Anaconda](https://www.continuum.io/downloads)

Environment installation:

```bash
	conda create -y --name prophyle \
		-c etetoolkit -c bioconda \
		python==3.4 ete3 ete3_external_apps bitarray \
		cmake parallel blast samtools
```

Environment activation:

```bash
	source activate prophyle
```

### Compile all programs

```bash
  make -C src
  make -C ext
```

### Download genomic libraries and simulate reads
```bash
  make -C library
```

### Download simulated reads
```bash
  make -C reads
```

### Custom taxonomic trees

Use [`bin/build_taxonomic_tree.py`](bin/build_taxonomic_tree.py) to build custom taxonomic trees starting from your database's fasta indexes and taxonomy files ([`library/Taxonomy`](library/Taxonomy) for more information). Taxonomic identifiers are assigned to the sequences first, and then the tree is built using [ETE Toolkit](http://etetoolkit.org/) and saved as newick format 1. Necessary node attributes are:

 * `name`: unique node name (format n[0-9]\*)
 * `taxid`: unique taxonomic identifier
 * `seqname`: names of the sequences sharing the same taxid, separated by @
 * `fastapath`: paths of the sequences' fasta files, separated by @ (absolute or relative from the main directory of the repository)
 * `infasta_offset`: positions where each sequence starts inside the corresponding fasta files, separated by @
 * `base_len`: length of each sequence, separated by @

Other optional attributes are `sci_name`, `named_lineage`, `lineage`, `rank` (more info [here](http://etetoolkit.org/docs/latest/tutorial/tutorial_ncbitaxonomy.html#automatic-tree-annotation-using-ncbi-taxonomy)).

## How to use

### Build the BWT-indexes

#### All experiments:

```bash
  make -C experiments -j 10
```

#### First experiment only

```bash
  make -C experiments/01* -j 10
```

