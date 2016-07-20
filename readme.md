# k-mer propagation experiments

[![Build Status](https://travis-ci.com/karel-brinda/MetaNG.svg?token=LzzDiQkWWqF4hBjZahmQ&branch=master)](https://travis-ci.com/karel-brinda/MetaNG)

## Getting started

### Prerequisities

* CMake 2.6+
* GCC 4.8+
* Boost 1.46+ 
* ZLib
* Python 3 with ete3 library
* SamTools

#### Recommended way of installation using [Anaconda](https://www.continuum.io/downloads)

Environment installation:

```bash
	conda create -y --name metang \
		-c etetoolkit -c bioconda \
		python==3.4 ete3 ete3_external_apps bitarray \
		cmake parallel blast
```

Environment activation:

```bash
	source activate metang
```

### Compile all programs

```bash
  make -C src
```

### Download genomic libraries and simulate reads
```bash
  make -C library
```

### Download simulated reads
```bash
  make -C reads
```

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
