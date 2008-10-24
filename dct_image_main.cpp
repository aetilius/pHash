#include <cstdio>
#include "CImg.h"
#include "pHash.h"

int main(int argc, char **argv){
    cimg_usage("pHash robust image hash program using Discrete Cosine Transform");
    const char *file1 = cimg_option("-f1",(char *)NULL,"name of first file");
    const char *file2 = cimg_option("-f2",(char *)NULL,"name of second file");
    const int thresh = cimg_option("-t",22,"threhold value");
    const char *msg = ph_about();
    puts(msg);

    printf("file: %s\n",file1);
    printf("file: %s\n",file2);
    printf("threshold is %d\n",thresh);
    int size = sizeof(ulong64);
    printf("hash has %d bytes\n",size);

    ulong64 hash1;
    if (ph_dct_imagehash(file1,hash1) < 0){
	return -1;
    }
    ulong64 hash2;
    if (ph_dct_imagehash(file2,hash2) < 0){
	return -1;
    }
    printf("hash1 is %llu\n",hash1);
    printf("hash2 is %llu\n",hash2);
    int hd = ph_hamming_distance(hash1,hash2);
    printf("hamming distance is %d\n",hd);
    if (hd < 0){
	printf("unable to get hamming distance\n");
	return -1;
    }
    if (hd > thresh)
	printf("images are different\n");
    else
	printf("images are same\n");

    return EXIT_SUCCESS;
}
