#include <cstdlib>
#include <cstdio>
#include <cassert>
#include "pHash.h"

using namespace std;

int main(int argc, char **argv){
	assert(argc == 3);
	const char *image_file1 = argv[1];
	const int expected_n1 = atoi(argv[2]);
	
	assert(image_file1 != NULL);
	
	int n1 = 0;
	ulong64 *hash1 = ph_dct_videohash(image_file1, n1);
	assert(hash1 != NULL);
	assert(n1 == expected_n1);

	printf("file: %s no. frames %d\n", image_file1, n1);

	return 0;
}
