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

// #define TEST_TAXA

using namespace std;
using namespace kraken;

static const char * INDEX = "database.idx";
static const char * DB = "database.kdb";
static const uint64_t NT_MASK = uint64_t(3) << 62;

char* kmerInterpreter(uint64_t kmer, uint8_t k) {
	char* res = new char[k];
	kmer <<= (sizeof(kmer) * 8) - (k * 2);
	for(int i=0; i<k; i++) {
		switch ((kmer & NT_MASK) >> 62) {
			case 0:
				res[i] = 'A';
				break;
			case 1:
				res[i] = 'C';
				break;
			case 2:
				res[i] = 'G';
				break;
			case 3:
				res[i] = 'T';
				break;
		}
		kmer <<= 2;
	}
	return res;
}

void test_taxa(KrakenDB db, uint64_t kmer, uint32_t taxid,
				uint64_t val_len, ofstream err_file)
{
	uint64_t last_bin_key = 0;
	int64_t min_pos = 1;
	int64_t max_pos = 0;
	uint32_t taxid_e = 0;
	uint32_t* kmer_pos = db.kmer_query(kmer,&last_bin_key,
							&min_pos,&max_pos,true);
	memcpy(&taxid_e, kmer_pos, val_len);
	err_file << "exp:\t" << taxid_e << endl;
	err_file << "act:\t" << taxid << endl;
}

int main() {

	#ifdef TEST_TAXA
	ofstream err_file;
	err_file.open("trans_kra.log");
	#endif
	
	QuickFile idx_file(INDEX);
	KrakenDBIndex idx(idx_file.ptr());
	idx_file.load_file();
	
	QuickFile db_file(DB);
	KrakenDB db(db_file.ptr());
	db.set_index(&idx);
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

    uint64_t* offsets = idx.get_array();
	uint64_t nt = idx.indexed_nt();
	uint64_t n_bins = (1ull << (2 * nt)) + 1;
	uint64_t bin_len = 64*key_ct/n_bins;
	
	uint64_t mask = (1ull << key_bits) -1;
	
	char* bin = (char*) malloc(bin_len*pair_size);
	uint64_t* kmers = (uint64_t*) calloc(bin_len, sizeof(uint64_t));
	uint32_t* taxa = (uint32_t*) calloc(bin_len, sizeof(uint32_t));
	uint64_t* temp_kmers;
	uint32_t* temp_taxa;
	uint64_t start = 0, stop = 0, len = 0;

	uint64_t kmer_id = 1;
	
	for(uint64_t* temp_offsets = offsets + 1;
			temp_offsets < offsets + n_bins;
			temp_offsets++)
	{
		
		temp_kmers = kmers;
		temp_taxa = taxa;
			
		start = stop;
		stop = *temp_offsets;
		len = stop - start;
		
		if (start < stop) {
			
			if (len > bin_len) {
				bin_len = 2*len;
				bin = (char*) realloc(bin, bin_len*pair_size);
				kmers = (uint64_t*) realloc(kmers, bin_len*sizeof(uint64_t));
				taxa = (uint32_t*) realloc(taxa, bin_len*sizeof(uint32_t));
				temp_kmers = kmers;
				temp_taxa = taxa;
			}
			
			memcpy(bin, pairs + (start * pair_size), len * pair_size);
			
			for(char* next_pair = bin;
				next_pair < bin + (len*pair_size);
				next_pair += pair_size)
			{
				
				memcpy(temp_kmers, next_pair, key_len);
				memcpy(temp_taxa, next_pair+key_len, val_len);
				*temp_kmers &= mask;
				
				#ifdef TEST_TAXA
				test_taxa(db, *temp_kmers, *temp_taxa, val_len, err_file);
				#endif
				
				temp_kmers++;
				temp_taxa++;
			}

			for(temp_kmers = kmers, temp_taxa = taxa;
				temp_kmers < kmers + len;
				temp_kmers++, temp_taxa++)
			{
				cout << ">kmer|" << kmer_id++ << "|taxid|" << *temp_taxa << endl;
				cout << kmerInterpreter(*temp_kmers, key_bits/2) << endl << endl;
			}
		}
	}

	free(bin);
	free(kmers);
	free(taxa);
	
	return 0;
}
