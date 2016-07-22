#include "kraken_headers.hpp"
#include "krakendb.hpp"
#include "quickfile.hpp"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

using namespace std;
using namespace kraken;

static const char * KRAKEN_INDEX2_STRING = "KRAKIX2";

string Input_DB_filename, Output_DB_filename, Index_filename, Taxa_Map_filename, Taxa_List_filename;
uint64_t Key_len = 8;
uint64_t Val_len = 4;
uint64_t Pair_size = 12;
uint64_t Key_ct;
uint64_t Taxa_count;

static map<uint64_t, uint64_t> map_taxa_fn(char* data);
static char* bin_and_sort_data(char* data, map<uint64_t, uint64_t> taxa_map, KrakenDBIndex idx);
static void make_lca_index(char* data, string index_filename, map<uint64_t, uint64_t> taxa_map);
static void parse_command_line(int argc, char **argv);
static void usage(int exit_code=EX_USAGE);

int main(int argc, char **argv) {

	parse_command_line(argc, argv);

	QuickFile input_db_file(Input_DB_filename);
	KrakenDB input_db(input_db_file.ptr());
	input_db_file.load_file();
	Key_len = input_db.get_key_len();
	Val_len = input_db.get_val_len();
	Key_ct = input_db.get_key_ct();
	Pair_size = Key_len + Val_len;
	char* pairs = (char*) input_db.get_pair_ptr();

	uint64_t skip_len = input_db.header_size();
	char* header = (char*) malloc(skip_len);
	memcpy(header, input_db_file.ptr(), skip_len);

	char* data = (char*) malloc(Pair_size*Key_ct);
	memcpy(data, pairs, Key_ct * Pair_size);
	map<uint64_t, uint64_t> taxa_map = map_taxa_fn(data);
	Taxa_count = taxa_map.size();
	make_lca_index(data, Index_filename, taxa_map);
	QuickFile index_file(Index_filename);
	KrakenDBIndex db_index(index_file.ptr());

	char* sort_data = bin_and_sort_data(data, taxa_map, db_index);

	FILE * output_file = fopen (Output_DB_filename.c_str(), "wb");
	fwrite(header, skip_len, 1, output_file);
	fwrite(sort_data, Pair_size, Key_ct, output_file);
	
	free(header);
	free(sort_data);
	fclose(output_file);
	return 0;
}

static map<uint64_t, uint64_t> map_taxa_fn(char* data) {
	map<uint64_t, uint64_t> taxa_map;
	uint32_t taxid;
	uint64_t count = 0;
	pair<map<uint64_t, uint64_t>::iterator,bool> ret;
	FILE* taxa_list_file = fopen(Taxa_List_filename.c_str(), "w");
	for(char* next_pair = data;
			next_pair < data + (Key_ct * Pair_size);
			next_pair += Pair_size)
	{
		memcpy(&taxid, next_pair+Key_len, Val_len);
		ret = taxa_map.insert(pair<uint64_t, uint64_t>((uint64_t)taxid,count));
		if (ret.second == true) {
			fprintf(taxa_list_file, "%u\n", taxid);
			count++;
		}
	}
	fclose(taxa_list_file);
	FILE* taxa_map_file = fopen(Taxa_Map_filename.c_str(), "wb");
	fwrite(&count, sizeof(uint64_t), 1, taxa_map_file);
	for(map<uint64_t, uint64_t>::iterator it = taxa_map.begin();
			it != taxa_map.end(); it++)
	{
		fwrite(&(it->first), sizeof(uint64_t), 1, taxa_map_file);
		fwrite(&(it->second), sizeof(uint64_t), 1, taxa_map_file);
	}
	fclose(taxa_map_file);
	return taxa_map;
}

static char* bin_and_sort_data(char* data, map<uint64_t, uint64_t> taxa_map, KrakenDBIndex idx) {
	uint64_t* offsets = idx.get_array();
	uint64_t* pos = (uint64_t*) malloc(sizeof(uint64_t) * (Taxa_count+1));
	char* sort_data = (char*) malloc(Pair_size * Key_ct);
	uint64_t taxid = 0;
	uint64_t i = 0;
	char* pair = (char*) malloc(Pair_size);
	char* pair_pos = NULL;
	
	memcpy(pos, offsets, (Taxa_count+1)*sizeof(uint64_t));

	for (char* next_pair = data;
			next_pair < data + (Key_ct*Pair_size);
			next_pair += Pair_size)
	{
		memcpy(pair, next_pair, Pair_size);
		memcpy(&taxid, next_pair + Key_len, Val_len);
		i = taxa_map.find(taxid)->second;
		pair_pos = sort_data + Pair_size*(pos[i]++);
		memcpy(pair_pos, pair, Pair_size);
	}

	free(pos);
	free(data);
	free(pair);
	
	return sort_data;
}

void make_lca_index(char* data, string index_filename, map<uint64_t, uint64_t> taxa_map) {
	uint64_t entries = taxa_map.size();
	uint64_t *bin_counts = new uint64_t[entries];
	uint8_t nt = 0;
	uint64_t taxid = 0;
	uint64_t tax_idx = 0;
	for (char* next_pair = data + Key_len;
			next_pair < data + (Key_ct*Pair_size);
			next_pair += Pair_size)
	{
		taxid = 0;
		tax_idx = 0;
		memcpy(&taxid, next_pair, Val_len);
		tax_idx = taxa_map.find(taxid)->second;
		bin_counts[tax_idx]++;
	}
	uint64_t *bin_offsets = new uint64_t[ entries + 1 ];
	bin_offsets[0] = 0;
	for (uint64_t i = 1; i <= entries; i++)
		bin_offsets[i] = bin_offsets[i-1] + bin_counts[i-1];

	QuickFile idx_file(index_filename, "w",
	strlen(KRAKEN_INDEX2_STRING) + 1 + sizeof(*bin_offsets) * (entries + 1));
	char *idx_ptr = idx_file.ptr();
	memcpy(idx_ptr, KRAKEN_INDEX2_STRING, strlen(KRAKEN_INDEX2_STRING));
	idx_ptr += strlen(KRAKEN_INDEX2_STRING);
	memcpy(idx_ptr++, &nt, 1);
	memcpy(idx_ptr, bin_offsets, sizeof(*bin_offsets) * (entries + 1));
}

void parse_command_line(int argc, char **argv) {
	int opt;

	if (argc > 1 && strcmp(argv[1], "-h") == 0)
		usage(0);
	while ((opt = getopt(argc, argv, "d:o:i:m:l:")) != -1) {
		switch (opt) {
			case 'd' :
				Input_DB_filename = optarg;
				break;
			case 'o' :
				Output_DB_filename = optarg;
				break;
			case 'i' :
				Index_filename = optarg;
				break;
			case 'm':
				Taxa_Map_filename = optarg;
				break;
			case 'l':
				Taxa_List_filename = optarg;
				break;
			default:
				usage();
				break;
		}
	}

	if (Input_DB_filename.empty() || Output_DB_filename.empty()
				|| Index_filename.empty() || Taxa_Map_filename.empty()
				|| Taxa_List_filename.empty())
		usage();
}

void usage(int exit_code) {
	cerr << "Usage: taxid_sorter -d input-db -o output-db -i output-idx -m taxa-map -l taxa-list\n";
	exit(exit_code);
}
