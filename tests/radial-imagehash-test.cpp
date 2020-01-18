#include <iostream>
#include <cassert>
#include <limits>
#include <cmath>
#include "pHash.h"

using namespace std;

int main(int argc, char **argv){
	assert(argc == 4);
	const char* image_file1 = argv[1];
	const char* image_file2 = argv[2];
	const double expected = atof(argv[3]);
	const int n_angles = 64;
	const double sigma = 1.0;
	const double gamma = 1.0;
	
	assert(image_file1 != NULL);
	assert(image_file2 != NULL);

	Digest digest1;
	int rc = ph_image_digest(image_file1, sigma, gamma, digest1, n_angles);
	assert(rc == 0);
	
	Digest digest2;
	rc = ph_image_digest(image_file2, sigma, gamma, digest2, n_angles);
	assert(rc == 0);

	double pcc;
	rc = ph_peakcrosscorr(digest1, digest2, pcc);
	assert(rc == 0);
	
	cout << "file: " << image_file1 << endl;
	cout << "file: " << image_file2 << endl;
	cout << "pcc = " << pcc << endl;
	
	double diff = fabs(pcc - expected);
	assert(diff <= 10*numeric_limits<float>::epsilon());
	
	return 0;
}
