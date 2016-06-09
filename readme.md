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
		snakemake samtools cmake parallel
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

Approx. time:
```
real    14m31.329s
user    3m20.748s
sys     1m40.644s
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

Approx. time:
```bash
real    21m4.152s
user    48m2.228s
sys     1m27.988s
```
