# k-mer propagation experiments

[![Build Status](https://travis-ci.com/karel-brinda/MetaNG.svg?token=LzzDiQkWWqF4hBjZahmQ&branch=master)](https://travis-ci.com/karel-brinda/MetaNG)

## Getting started

### Prerequisities

* CMake 2.6+
* GCC 4.8+
* Boost 1.46+ 
* ZLib
* Python 3 with ete3 library
* GNU time (on OS X using ``brew install gnu-time``)

#### Recommended way of installation using [Anaconda](https://www.continuum.io/downloads)

Environment installation:

```bash
	conda create -y --name metang \
		-c etetoolkit -c bioconda \
		python==3.4 ete3 ete3_external_apps \
		cmake parallel
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

Approx. time:
```
real    14m31.329s
user    3m20.748s
sys     1m40.644s
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

Approx. time:
```bash
real    21m4.152s
user    48m2.228s
sys     1m27.988s
```

#### Redefinition of ``time``

Useful when the default ``time`` installed in the system does not support verbose mode ``-v``.

Example:

```bash
make -j 7 TIME=~/.linuxbrew/bin/time -C experiments/
```
