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

using std::string;
using std::vector;
using namespace std;
using namespace kraken;

static const char * INDEX = "database.idx";
static const char * DB = "database.kdb";
static const uint64_t MASK = uint64_t(3) << 62;

char* kmerInterpreter(uint64_t kmer, uint8_t k) {
	char* res = new char[k];
	kmer <<= sizeof(kmer) * 8 - (k * 2);
	for(int i=0; i<k; i++) {
		switch ((kmer & MASK) >> 62) {
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

int main() {
	ofstream out_file;
	out_file.open("out");
	QuickFile idx_file(INDEX);
	KrakenDBIndex idx(idx_file.ptr());
	idx_file.load_file();
	QuickFile db_file(DB);
	KrakenDB db(db_file.ptr());
	db.set_index(&idx);
	// uint64_t key_len = db.get_key_len();
	// uint64_t val_len = db.get_val_len();
	// uint64_t key_ct = db.get_key_ct();
	// uint64_t skip_len = db.header_size();
	uint64_t* bins = idx.get_array();
	uint64_t nt = idx.indexed_nt();
	uint64_t sz = sizeof(nt);
	for(uint64_t i = 0; i < 1 + (1ull << (nt * 2)); i++) {
		out_file << bins[i] << endl;
	}
	
	return 0;
}
