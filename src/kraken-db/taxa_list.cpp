#include "kraken_headers.hpp"

using namespace std;

static void parse_command_line(int argc, char **argv);
static void usage(int exit_code=EX_USAGE);

string Taxa_Map_filename, Output_filename;

int main(int argc, char **argv) {

	parse_command_line(argc, argv);
	
	FILE* taxa_map_file = fopen(Taxa_Map_filename.c_str(), "rb");
	uint64_t taxa_count;
	fread(&taxa_count, sizeof(uint64_t), 1, taxa_map_file);
	uint64_t* taxa_map = (uint64_t*) malloc(2*taxa_count*sizeof(uint64_t));
	fread(taxa_map, sizeof(uint64_t), 2*taxa_count, taxa_map_file);
	fclose(taxa_map_file);

	FILE* output_file = fopen(Output_filename.c_str(), "w");
	for (uint64_t* next_entry = taxa_map;
			next_entry < taxa_map + 2*taxa_count;
			next_entry += 2)
	{
		fprintf(output_file, "%lu\n", *next_entry);
	}
	fclose(output_file);
	free(taxa_map);

	return 0;
}

void parse_command_line(int argc, char **argv) {
	int opt;
	
	if (argc > 1 && strcmp(argv[1], "-h") == 0)
		usage(0);
	while ((opt = getopt(argc, argv, "m:o:")) != -1) {
		switch (opt) {
			case 'm':
				Taxa_Map_filename = optarg;
				break;
			case 'o':
				Output_filename = optarg;
				break;
			default:
				usage();
				break;
		}
	}
	if (Taxa_Map_filename.empty() || Output_filename.empty())
		usage();
}

void usage(int exit_code) {
	cerr << "Usage: taxa_list -m taxa-map -o output-file\n";
	exit(exit_code);
}
