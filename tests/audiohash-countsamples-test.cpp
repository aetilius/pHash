#include <cstdlib>
#include <iostream>
#include <cassert>
#include <audiophash.h>

using namespace std;

int main(int argc, char **argv){
	assert(argc == 5);
	const char* image_file1 = argv[1];
	const char* image_file2 = argv[2];
	const int n_expected1 = atoi(argv[3]);
	const int n_expected2 = atoi(argv[4]);
	const int sr = 8000;
	const int n_channels = 1;
	
	assert(image_file1 != NULL);
	assert(image_file2 != NULL);

	int n1 = ph_count_samples(image_file1, sr, n_channels);
	assert(n1 == n_expected1);

	int n2 = ph_count_samples(image_file2, sr, n_channels);
	assert(n2 == n_expected2);
	
	cout << "file: " << image_file1 << " samples " << n1 << endl;
	cout << "file: " << image_file2 << " samples " << n2 << endl;
	
	return 0;
}
