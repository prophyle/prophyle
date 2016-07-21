#include "../kraken_headers.hpp"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

using namespace std;

static const char * DATABASE_FILE_TYPE = "JFLISTDN";

static void parse_command_line(int argc, char **argv);
static void usage(int exit_code=EX_USAGE);

string Input_DB_filename, Output_DB_filename;
uint64_t Kmers_Count = 0;
uint64_t Taxa_Count = 0;

int main (int argc, char **argv) {

	parse_command_line(argc, argv);

	uint64_t key_bits = 62;
	uint64_t key_len = 8;
	uint64_t val_len = 4;
	uint64_t pair_size = key_len + val_len;
	uint64_t skip_len = 72 + 2 * (4 + 8 * key_bits);
	char* header = (char*) calloc(1, skip_len);
	memcpy(header, DATABASE_FILE_TYPE, 8);
	memcpy(header + 8, &key_bits, 8);
	memcpy(header + 16, &val_len, 8);
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
	while ((opt = getopt(argc, argv, "k:t:o:")) != -1) {
		switch (opt) {
			case 'k':
				sscanf(optarg, "%lu", &Kmers_Count);
				break;
			case 't':
				sscanf(optarg, "%lu", &Taxa_Count);
				break;
			case 'o':
				Output_DB_filename = optarg;
				break;
			default:
				usage();
				break;
		}
	}
	if (Kmers_Count == 0 || Taxa_Count == 0 || Output_DB_filename.empty())
		usage();
}

void usage(int exit_code) {
	cerr << "Usage: build_test_db -k kmers-count -t taxa-count -o output-db\n";
	exit(exit_code);
}
