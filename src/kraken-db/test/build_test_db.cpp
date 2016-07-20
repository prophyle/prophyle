#include "kraken_headers.hpp"
#include "quickfile.hpp"
#include "krakendb.hpp"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

using namespace std;
using namespace kraken;

static void parse_command_line(int argc, char **argv);
static void usage(int exit_code=EX_USAGE);

string Input_DB_filename, Output_DB_filename;
uint64_t Kmers_Count = 0;
uint64_t Taxa_Count = 0;

int main (int argc, char **argv) {

	parse_command_line(argc, argv);

	QuickFile input_db_file(Input_DB_filename);
	KrakenDB input_db(input_db_file.ptr());
	uint64_t key_len = input_db.get_key_len();
	uint64_t key_bits = input_db.get_key_bits();
	uint64_t val_len = input_db.get_val_len();
	uint64_t pair_size = key_len + val_len;
	uint64_t k_mask = (1ull << key_bits) -1;

	uint64_t skip_len = input_db.header_size();
	char* header = (char*) malloc(skip_len);
	memcpy(header, input_db_file.ptr(), skip_len);
	memcpy(header + 48, &Kmers_Count, 8);

	uint64_t kmer = 0;
	uint32_t taxid = 0;
	char* data = (char*) malloc(pair_size*Kmers_Count);
	char* pair_pos = data;
    
	for(uint64_t i = 0; i < Kmers_Count; i++) {
		taxid = i%Taxa_Count + 1;
		kmer = i;
		memcpy(pair_pos, &kmer, key_len);
		memcpy(pair_pos + key_len, &taxid, val_len);
		pair_pos += pair_size;
	}

	FILE * output_file = fopen (Output_DB_filename.c_str(), "w");
	fwrite(header, skip_len, 1, output_file);
	fwrite(data, pair_size, Kmers_Count, output_file);

	fclose(output_file);
	free(header);
	free(data);
}

void parse_command_line(int argc, char **argv) {
	int opt;
	
	if (argc > 1 && strcmp(argv[1], "-h") == 0)
		usage(0);
	while ((opt = getopt(argc, argv, "k:t:d:o:")) != -1) {
		switch (opt) {
			case 'k':
				sscanf(optarg, "%lu", &Kmers_Count);
				break;
			case 't':
				sscanf(optarg, "%lu", &Taxa_Count);
				break;
			case 'd':
				Input_DB_filename = optarg;
				break;
			case 'o':
				Output_DB_filename = optarg;
				break;
			default:
				usage();
				break;
		}
	}
	if (Kmers_Count == 0 || Taxa_Count == 0 ||
			Input_DB_filename.empty() || Output_DB_filename.empty())
		usage();
}

void usage(int exit_code) {
	cerr << "Usage: build_test_db -k kmers-count -t taxa-count -d input-db -o output-db\n";
	exit(exit_code);
}
