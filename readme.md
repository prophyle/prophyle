# k-mer propagation experiments

## Prerequisities

* CMake
* BWA
* C++ with Boost
* Python 3

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
