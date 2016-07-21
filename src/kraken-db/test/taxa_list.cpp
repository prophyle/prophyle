#include "kraken_headers.hpp"

using namespace std;

int main(int argc, char **argv) {

	const char* taxa_map_filename = argv[1];
	FILE* taxa_map_file = fopen(taxa_map_filename, "rb");
	uint64_t taxa_count;
	fread(&taxa_count, sizeof(uint64_t), 1, taxa_map_file);
	uint64_t* taxa_map = (uint64_t*) malloc(2*taxa_count*sizeof(uint64_t));
	fread(taxa_map, sizeof(uint64_t), 2*taxa_count, taxa_map_file);
	fclose(taxa_map_file);

	for (uint64_t* next_entry = taxa_map;
			next_entry < taxa_map + 2*taxa_count;
			next_entry += 2)
	{
		printf("%lu\n", *next_entry);
	}
	free(taxa_map);

	return 0;
}
