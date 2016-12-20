## Paper outline

Target journal: PLOS Comp Bio?? Genome Reseach (as Centrifuge)?? Nature Scientific Reports?? Nature Communications (as Kaiju)?? Nature Methods?? Nature Biotechnology (as Sequence Bloom Trees)?? BMC Genome Biology (as Kraken)????

Proposed title: ProPhyle: accurate and resource-frugal phylogeny-based metagenomic classification

Goals of paper: 
* announce ProPhyle
* demonstrate convincingly that (a) more accurate, and (b) much more memory efficient, than Kraken (and Centrifuge)

Model articles:
* ??

## Plan (structure taken from PLOS)

#### Abstract
##### Keywords: next-generation sequencing, metagenomic classification, phylogeny, alignment-free sequence comparison

#### Author summary

### Introduction

### Results

### Discussion

Explain carefully how ProPhyle relates to Kraken and Centrifuge.

Kraken: 
* ProPhyle index is more expressive, as it characterizes precisely in which genomes each k-mer occurs
* as a consequence, ProPhyle tends to classify lower in the tree, which should lead to a better sensitivity
* yet our index is small, much smaller than the one of Kraken
* Kraken's measures (given in the paper) are misleading because they mix all genomes together. We expect to be uniformely better for individual genomes
* what about the assignment algorithm?

Centrifuge:
* ProPhyle shares with Centrifuge similar ideas: "propagating" common sequences, using BWT, outputting several assignments for a read (however, in their experiments, only uniquely assigned reads are considered)
* However, Centrifuge makes various empirical steps: identifying shared sequences (done with the external Nucmer tool), neighbour-joining merging algorithm, scoring species, reducing the number of assignments by their top-down propagation. All these steps are very empirical or even naive. 
* In contrast, ProPhyle works at the level of k-mers, in a "lossless" fashion
* Note that, in contrast to Kraken, Centrifuge does not really maps reads to the tree, but classifies them to (one or several) genomes and then possibly propagates bottom-up to some internal node. ProPhyle (and Kraken) directly map every read to a (possibly internal) tree node
* Centrifuge index takes 4.2Gb for ~4300 bacterial genomes. **Can we compete?**
* *question: what BWT-index implementation is used in Centrifuge?*

Compared to both Kraken and Centrifuge:
* our program is "self-contained" i.e. does not use external tools (e.g. Centrifuge uses Nucmer for constructing the index)
* since we work with Newick tree format, we can use any phylogeny, not necessarily NCBI
* [*other reasons why ProPhyle has been programmed more professionally*]
* what about the classification speed (?)

### Methods and Models

### Acknowledgements

### References
