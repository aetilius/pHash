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

	BMBHash hash1, hash2;
	int rc = ph_bmb_imagehash(image_file1, hash1);
	assert(rc == 0);
	
	rc = ph_bmb_imagehash(image_file2, hash2);
	assert(rc == 0);
	
	int d = ph_bmb_distance(hash1, hash2);
	assert(d >= 0);

	ph_bmb_free(hash1);
	ph_bmb_free(hash2);
	
	printf("file: %s hash[%d]\n", image_file1, hash1.bytelength);
	printf("file: %s hash[%d]\n", image_file2, hash2.bytelength);
	printf("d = %d\n", d);
	
	return 0;
}
