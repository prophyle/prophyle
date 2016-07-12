#include "kraken_headers.hpp"
#include "krakendb.hpp"
#include "quickfile.hpp"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

using namespace kraken;

static const char * DB = "database.kdb";
static const uint64_t NT_MASK = uint64_t(3) << 62;

char* kmerInterpreter(uint64_t kmer, uint8_t k, char* trans_kmer) {
	kmer <<= (sizeof(kmer) * 8) - (k * 2);
	for(int i=0; i<k; i++) {
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

int main() {

	QuickFile db_file(DB);
	KrakenDB db(db_file.ptr());
	db_file.load_file();
	
	// Return pointer to start of pairs
	char* pairs = db.get_pair_ptr();
	// how many bits are in each key?
    uint64_t key_bits = db.get_key_bits();
    // how many bytes does each key occupy?
    uint64_t key_len = db.get_key_len();
    // how many bytes does each value occupy?
    uint64_t val_len = db.get_val_len();
    // how many key/value pairs are there?
    uint64_t key_ct = db.get_key_ct();
    // how many bytes does each pair occupy?
    uint64_t pair_size = db.pair_size();

    uint64_t mask = (1ull << key_bits) -1;
	
	char* bin = (char*) malloc(key_ct * pair_size);
	char *trans_kmer = (char*) malloc(key_bits/2);
	uint64_t kmer;
	uint32_t taxid;
	unsigned long long int kmer_id = 1;
	
	memcpy(bin, pairs, key_ct * pair_size);

	for(char* next_pair = pairs;
			next_pair < pairs + (key_ct * pair_size);
			next_pair += pair_size)
	{
		memcpy(&kmer, next_pair, key_len);
		memcpy(&taxid, next_pair+key_len, val_len);
		kmer &= mask;
		kmerInterpreter(kmer, key_bits/2, trans_kmer);
		printf(">kmer|%llu|taxid|%lu\n%s\n\n", kmer_id++, (unsigned long)taxid, trans_kmer);
	}

	free(bin);
	free(trans_kmer);	
	return 0;
}
