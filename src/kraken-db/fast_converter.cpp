#include "kraken_headers.hpp"
#include "krakendb.hpp"
#include "quickfile.hpp"

#include <sys/stat.h>
#include <sys/mman.h> 
#include <errno.h>
#include <string.h>
#include <stdarg.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

using namespace std;
using namespace kraken;

char* kmer_interpreter(uint64_t kmer, uint64_t k, char* trans_kmer);
static void parse_command_line(int argc, char **argv);
static void usage(int exit_code=EX_USAGE);

string Input_DB_filename, Output_filename;
static const uint64_t NT_MASK = uint64_t(3) << 62;

int main(int argc, char **argv) {

	parse_command_line(argc, argv);

	QuickFile db_file(Input_DB_filename);
	KrakenDB db(db_file.ptr());
	db_file.load_file();
	
	char* pairs = db.get_pair_ptr();
	uint64_t key_bits = db.get_key_bits();
    uint64_t key_len = db.get_key_len();
    uint64_t val_len = db.get_val_len();
    uint64_t key_ct = db.get_key_ct();
    uint64_t pair_size = db.pair_size();

    uint64_t mask = (1ull << key_bits) -1;
	
	char* bin = (char*) malloc(key_ct * pair_size);
	char *trans_kmer = (char*) malloc(key_bits/2);
	uint64_t kmer;
	uint32_t taxid;
	unsigned long long int kmer_id = 1;
	
	memcpy(bin, pairs, key_ct * pair_size);
	FILE * output_file = fopen(Output_filename.c_str(), "w");
	
	for(char* next_pair = pairs;
			next_pair < pairs + (key_ct * pair_size);
			next_pair += pair_size)
	{
		memcpy(&kmer, next_pair, key_len);
		memcpy(&taxid, next_pair+key_len, val_len);
		kmer &= mask;
		kmer_interpreter(kmer, key_bits, trans_kmer);
		fprintf(output_file, ">kmer|%llu|taxid|%lu\n%s\n\n", kmer_id++, (unsigned long)taxid, trans_kmer);
	}

	fclose(output_file);
	free(bin);
	free(trans_kmer);

	return 0;
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
	while ((opt = getopt(argc, argv, "d:o:")) != -1) {
		switch (opt) {
			case 'd':
				Input_DB_filename = optarg;
				break;
			case 'o':
				Output_filename = optarg;
				break;
			default:
				usage();
				break;
		}
	}
	if (Input_DB_filename.empty() || Output_filename.empty())
		usage();
}

void usage(int exit_code) {
	cerr << "Usage: fast_converter -d input-db -o output-file\n";
	exit(exit_code);
}
