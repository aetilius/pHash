#include <cstdlib>
#include <iostream>
#include <cmath>
#include <cassert>
#include <audiophash.h>

using namespace std;

int main(int argc, char **argv){
	assert(argc == 8);
	const char* image_file1 = argv[1];
	const char* image_file2 = argv[2];
	const int expected_nframes1 = atoi(argv[3]);
	const int expected_nframes2 = atoi(argv[4]);
	const float expected_max_cs = atof(argv[5]);
	const int expected_pos = atoi(argv[6]);
	const int expected_Nc = atoi(argv[7]);
	const int sr = 8000;
	const int n_channels = 1;
	const float threshold = 0.25;
	const int bs = 16;
	
	assert(image_file1 != NULL);
	assert(image_file2 != NULL);

	int n1;
	float *buf1 = ph_readaudio(image_file1, sr, n_channels, NULL, n1, 0);
	assert(buf1 != NULL);
	assert(n1 >= 0);
	
	int n2;
	float *buf2 = ph_readaudio(image_file2, sr, n_channels, NULL, n2, 0);
	assert(buf2 != NULL);
    assert(n2 >= 0);
	
	int n_frames1;
	uint32_t *hash1 = ph_audiohash(buf1, n1, sr, n_frames1);
	assert(hash1 != NULL);
	assert(n_frames1 == expected_nframes1);
	
	int n_frames2;
	uint32_t *hash2 = ph_audiohash(buf2, n2, sr, n_frames2);
	assert(hash2 != NULL);
	assert(n_frames2 == expected_nframes2);

	int Nc;
	double *cs = ph_audio_distance_ber(hash1, n_frames1, hash2, n_frames2, threshold, bs, Nc);
	assert(cs != NULL);
	assert(Nc == expected_Nc);

	double max_cs = 0;
	int pos = -1;
	for (int i=0;i < Nc;i++){
		if (cs[i] > max_cs){
			max_cs = cs[i];
			pos = i;
		}
	}
	assert(pos == expected_pos);
	assert(fabs(max_cs - expected_max_cs) <= 100*numeric_limits<float>::epsilon());
	
	cout << "file: " << image_file1 << " frames " << n_frames1 << endl;
	cout << "file: " << image_file2 << " frames " << n_frames2 << endl;
	cout << "confidence vector " << Nc << endl;
	cout << "confidence peak " << max_cs << " at position " << pos << endl;
	
	free(buf1);
	free(buf2);
	free(hash1);
	free(hash2);
	
	return 0;
}
