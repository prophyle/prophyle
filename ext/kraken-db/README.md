## Krakend-DB tools

Set of tools to compare MetaNG results to kraken's.

### `taxid_sorter`

Sorts a kraken-like database by taxid (instead of minimizer) and creates an index for it.

Options:

* `-d`: input database file. If k-mers in each bin need to be lexicographically sorted, use kraken's db_sort script first;
* `-o`: output file for the sorted database;
* `-i`: output file used to store the new index, which will contain the offsets of the taxids' "bins";
* `-m`: output file for the map containg (taxid, bin) associations, where bin is the index of the taxid in the offsets array;
* `-l`: output file for the list of tax 

### `get_kmers_by_taxid`

Extracts all the k-mers associated to the given taxid and prints them to stdout.

Options:

* `-t`: taxonomic identifier whose k-mers should be extracted;
* `-d`: input database file, sorted by `taxid_sorter`;
* `-i`: input index file, produced by `taxid_sorter`;
* `-m`: taxid - bin map, produced by `taxid_sorter`.

### `kraken_assembler.sh`

Calls MetaNG's assembler program using the ouput of `get_kmers_by_taxid` as argument, creating contigs for each taxid. It assumes there is a directory named `kraken` in the script's one, and that it contains `get_kmers_by_taxid`.

Options:

* `-l`: list of taxonomic identifiers whose k-mers should be extracted and assembled (output of `taxid_sorter -l`);
* `-m`: taxid - bin map, produced by `taxid_sorter`;
* `-d`: input database file, sorted by `taxid_sorter`;
* `-i`: input index file, produced by `taxid_sorter`;
* `-f`: directory in which the MetaNG's index will be stored;
* `-k`: k-mers lenght.
