#include <stdint.h>
#include <stdio.h>
#include "math.h"
#include "pHash.h"

using namespace cimg_library;

int main(int argc, char **argv) {
    if (argc < 3) {
        printf(" not enough input args\n");
        exit(1);
    }

    const char *img1 = argv[1];
    const char *img2 = argv[2];

    int alpha = 2;
    int level = 1;

    int hashlen1, hashlen2;

    printf("file1: %s\n", img1);
    uint8_t *hash1 = ph_mh_imagehash(img1, hashlen1, alpha, level);
    printf("file2: %s\n", img2);
    uint8_t *hash2 = ph_mh_imagehash(img2, hashlen2, alpha, level);
    if (hash2 == NULL || hash1 == NULL) {
        printf("Unable to generate hash");
        return (1);
    }
    double dist = ph_hammingdistance2(hash1, hashlen1, hash2, hashlen2);
    printf("distance = %f\n", dist);
    printf("-------------\n");
    printf("done\n");

    free(hash1);
    free(hash2);
    return 0;
}
