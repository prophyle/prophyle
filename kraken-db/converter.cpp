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

static void check (int test, const char * message, ...)
{
    if (test) {
        va_list args;
        va_start (args, message);
        vfprintf (stderr, message, args);
        va_end (args);
        fprintf (stderr, "\n");
        exit (EXIT_FAILURE);
    }
}

size_t getFilesize(const char* filename) {
	struct stat st;
	if(stat(filename, &st) != 0) {
		return 0;
	}
	return st.st_size;   
}

void openFileRO(const char* filename, const char **fptr) {
	int fd = open(filename, O_RDONLY);
	check (fd < 0, "open %s failed: %s", filename, strerror (errno));
	size_t size = getFilesize(filename);
	check (size == 0, "stat %s failed: %s", filename, strerror (errno));
	*fptr = (char*) mmap(0, size, PROT_READ, MAP_PRIVATE, fd, 0);
	check (*fptr == MAP_FAILED, "mmap %s failed: %s",
           filename, strerror (errno));
	//if (fptr == MAP_FAILED)
	//	err(EX_OSERR, "unable to mmap %s", filename);
}

int main() {
	const char* idx_file;
	openFileRO(INDEX, &idx_file);
	char* kraken_idx_file = (char *) idx_file;
	KrakenDBIndex idx(kraken_idx_file);
}
