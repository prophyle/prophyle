# k-mer propagation experiments

## Getting started

### Prerequisities

* GIT
* CMake
* C++ with Boost
* Python 3 with ete3 library, Snakemake, and RNFtools
* SamTools

#### Recommended way of installation using [Anaconda](https://www.continuum.io/downloads)

Environment installation:

```bash
	conda create -y --name metang \
		-c etetoolkit -c bioconda \
		python==3.4 ete3 ete3_external_apps \
		snakemake samtools git cmake parallel
```

Environment activation:

```bash
	source activate metang
```

RNFtools installation (in the activated environment)

```bash
	pip install git+git://github.com/karel-brinda/rnftools
```

### Compile all programs

```bash
  make -C src
```

### Download genomic libraries and simulate reads
```bash
  make -C libraries
```

## How to use

### Build the BWT-indexes

For all experiments:

```bash
  make -C experiments -j 10
```

For first experiment only (quick testing):

```bash
  make -C experiments/01* -j 10
```
