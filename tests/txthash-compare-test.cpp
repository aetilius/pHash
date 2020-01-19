#include <cstdlib>
#include <cstdio>
#include <cassert>
#include "pHash.h"

using namespace std;

int main(int argc, char **argv){
	assert(argc == 6);
	const char* text_file1 = argv[1];
	const char* text_file2 = argv[2];
	const int expected_n1 = atoi(argv[3]);
	const int expected_n2 = atoi(argv[4]);
	const int expected_matches = atoi(argv[5]);
	
	assert(text_file1 != NULL);
	assert(text_file2 != NULL);

	int n1 = 0;
	TxtHashPoint *txthash1 = ph_texthash(text_file1, &n1);
	assert(txthash1 != NULL);
	assert(n1 == expected_n1);

	int n2 = 0;
	TxtHashPoint *txthash2 = ph_texthash(text_file2, &n2);
	assert(txthash2 != NULL);
	assert(n2 == expected_n2);

	int n_matches = 0;
	TxtMatch *matches = ph_compare_text_hashes(txthash1, n1, txthash2, n2, &n_matches);
	assert(matches != NULL);
	assert(n_matches == expected_matches);

	printf("file: %s no. kgrams %d\n", text_file1, n1);
	printf("file: %s no. kgrams %d\n", text_file2, n2);
	printf("no. matches %d\n", n_matches);
	for (int i=0;i<n_matches;i++){
		printf("(%d) first index %ld second index %d  length %d\n",
			   i, matches[i].first_index, matches[i].second_index, matches[i].length);
	}
	
	return 0;
}
