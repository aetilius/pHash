#include <iostream>
#include <cassert>
#include "pHash.h"

using namespace std;

int main(int argc, char **argv){
	assert(argc == 4);
	const char* image_file1 = argv[1];
	const char* image_file2 = argv[2];
	const int expected = atoi(argv[3]);
	const float alpha = 2.0;
	const float lvl = 1.0;
	
	assert(image_file1 != NULL);
	assert(image_file2 != NULL);

	int N1;
	uint8_t *hash1 = ph_mh_imagehash(image_file1, N1, alpha, lvl);
	assert(hash1 != NULL);
	int N2;
	uint8_t *hash2 = ph_mh_imagehash(image_file2, N2, alpha, lvl);
	assert(hash2 != NULL);

	int d = ph_hammingdistance2(hash1, N1, hash2, N2);
	assert(d == expected);

	cout << "file: " << image_file1 << "hash[" << N1 << "]" << endl;
	cout << "file: " << image_file2 << "hash[" << N2 << "]" << endl;
	cout << "d = " << d << endl;
	
	return 0;
}
