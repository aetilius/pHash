#include <cstdio>
#include "pHash.h"
#include "cimgffmpeg.h"

using namespace cimg_library;


int main(int argc, char **argv){
    const char *file1 = cimg_option("-f1",(char *)NULL,"file name of first video");
    const char *file2 = cimg_option("-f2",(char *)NULL,"file name of second video");
    const char *msg = ph_about();
    puts(msg);
    printf("file1: %s\n",file1);
    printf("file2: %s\n",file2);

    ulong64 hash1;
    if (ph_dct_videohash(file1,hash1) < 0){
	printf("unable to calculate hash\n");
        exit(1);
    }

    ulong64 hash2;
    if (ph_dct_videohash(file2,hash2) < 0){
	printf("unable to calculate hash\n");
	exit(1);
    }

    int dist = ph_hamming_distance(hash1,hash2);
    printf("distance is %d\n",dist);

    return EXIT_SUCCESS;
}
