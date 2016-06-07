# k-mer propagation experiments

## Prerequisities

* GIT
* CMake
* C++ with Boost
* Python 3 with ete3 library, Snakemake, and RNFtools
* SamTools

### Recommended way of installation using [Anaconda](https://www.continuum.io/downloads)

Environment installation:

```bash
	conda create -y --name metang -c etetoolkit -c bioconda python==3.4 ete3 ete3_external_apps snakemake samtools git cmake
```

Environment activation:

```bash
	source activate metang
```

RNFtools installation (in the activated environment)

```bash
	pip install git+http://github.com/karel-brinda/rnftools
```

## Compilation of all programs

```bash
  make -C src
```

## Downloading genomic libraries
```bash
  make -C libraries
```

## Building the indexes

```bash
  make -C experiments -j 10
```

For a quick experiment:

```bash
  make -C experiments/01* -j 10
```
