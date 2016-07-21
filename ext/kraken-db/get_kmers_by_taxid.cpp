#include "kraken_headers.hpp"
#include "krakendb.hpp"
#include "quickfile.hpp"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

using namespace std;
using namespace kraken;

string Input_DB_filename, Index_filename, Taxa_Map_filename;
uint64_t Taxid = 0;
static const uint64_t NT_MASK = uint64_t(3) << 62;

static void convert_bin(KrakenDB db, KrakenDBIndex idx, uint64_t i);
uint64_t find_index(uint64_t* taxa_map, uint64_t taxa_count);
static char* kmer_interpreter(uint64_t kmer, uint64_t k, char* trans_kmer);
static void parse_command_line(int argc, char **argv);
static void usage(int exit_code=EX_USAGE);

int main(int argc, char **argv) {

	parse_command_line(argc, argv);
	
	QuickFile input_db_file(Input_DB_filename);
	KrakenDB input_db(input_db_file.ptr());
	QuickFile index_file(Index_filename);
	KrakenDBIndex db_index(index_file.ptr());
	
	FILE* taxa_map_file = fopen(Taxa_Map_filename.c_str(), "rb");
	uint64_t taxa_count;
	fread(&taxa_count, sizeof(uint64_t), 1, taxa_map_file);
	uint64_t* taxa_map = (uint64_t*) malloc(2*taxa_count*sizeof(uint64_t));
	fread(taxa_map, sizeof(uint64_t), 2*taxa_count, taxa_map_file);
	fclose(taxa_map_file);

	uint64_t bin_idx = find_index(taxa_map, taxa_count);
	free(taxa_map);

	convert_bin(input_db, db_index, bin_idx);

	return 0;
}

void convert_bin(KrakenDB db, KrakenDBIndex idx, uint64_t bin_idx) {
	uint64_t key_len = db.get_key_len();
	uint64_t key_bits = db.get_key_bits();
	uint64_t val_len = db.get_val_len();
	uint64_t key_ct = db.get_key_ct();
	char* pairs = (char*) db.get_pair_ptr();
	uint64_t pair_size = key_len + val_len;
	uint64_t* offsets = idx.get_array();
	uint64_t k_mask = (1ull << key_bits) -1;

	uint64_t start = offsets[bin_idx];
	uint64_t stop = offsets[bin_idx+1];
	uint64_t len = stop - start;

	char* bin = (char*) malloc(len*pair_size);
	memcpy(bin, pairs + start*pair_size, len*pair_size);
	
	uint64_t temp_kmer;
	uint32_t temp_taxid;
	uint64_t kmer_id = 0;
	char* trans_kmer = (char*) malloc(key_bits/2);
	
	for(char* next_pair = bin;
		next_pair < bin + (len*pair_size);
		next_pair += pair_size)
	{
		memcpy(&temp_kmer, next_pair, key_len);
		memcpy(&temp_taxid, next_pair+key_len, val_len);
		temp_kmer &= k_mask;
		kmer_interpreter(temp_kmer, key_bits, trans_kmer);
		printf(">kmer|%lu|taxid|%u\n%s\n",
			++kmer_id, temp_taxid, trans_kmer);
	}

	free(bin);
	free(trans_kmer);
}

uint64_t find_index(uint64_t* taxa_map, uint64_t taxa_count) {
	uint64_t i = 0;
	for (uint64_t* next_entry = taxa_map;
			next_entry < taxa_map + 2*taxa_count;
			next_entry += 2)
	{
		if (Taxid == *next_entry) {
			i = *(++next_entry);
			break;
		}
	}
	return i;
}

char* kmer_interpreter(uint64_t kmer, uint64_t k, char* trans_kmer) {
	kmer <<= (sizeof(kmer) * 8) - k;
	for(int i=0; i<k/2; i++) {
		switch ((kmer & NT_MASK) >> 62) {
			case 0:
				trans_kmer[i] = 'A';
				break;
			case 1:
				trans_kmer[i] = 'C';
				break;
			case 2:
				trans_kmer[i] = 'G';
				break;
			case 3:
				trans_kmer[i] = 'T';
				break;
		}
		kmer <<= 2;
	}
}

void parse_command_line(int argc, char **argv) {
	int opt;
	
	if (argc > 1 && strcmp(argv[1], "-h") == 0)
		usage(0);
	while ((opt = getopt(argc, argv, "t:d:i:m:")) != -1) {
		switch (opt) {
			case 't':
				sscanf(optarg, "%lu", &Taxid);
				break;
			case 'd':
				Input_DB_filename = optarg;
				break;
			case 'i':
				Index_filename = optarg;
				break;
			case 'm':
				Taxa_Map_filename = optarg;
				break;
			default:
				usage();
				break;
		}
	}
	if (Taxid == 0 || Input_DB_filename.empty() || Index_filename.empty() || Taxa_Map_filename.empty())
		usage();
}

void usage(int exit_code) {
	cerr << "Usage: get_kmers_by_taxid -t taxid -d input-db -i input-idx -m taxa-map\n";
	exit(exit_code);
}
