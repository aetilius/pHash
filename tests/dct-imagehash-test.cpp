#include <cstdlib>
#include <cstdio>
#include <cassert>
#include "pHash.h"

using namespace std;

int main(int argc, char **argv){
	assert(argc == 4);
	const char* image_file1 = argv[1];
	const char* image_file2 = argv[2];
	const int expected = atoi(argv[3]);
	
	assert(image_file1 != NULL);
	assert(image_file2 != NULL);

	ulong64 hash1;
	int rc = ph_dct_imagehash(image_file1, hash1);
	assert(rc == 0);

	ulong64 hash2;
	rc = ph_dct_imagehash(image_file2, hash2);
	assert(rc == 0);
	
	int d = ph_hamming_distance(hash1, hash2);
	assert(d == expected);

	printf("file: %x %s\n", hash1, image_file1);
	printf("file: %x %s\n", hash2, image_file2);
	printf("d = %d\n", d);
	
	return 0;
}
