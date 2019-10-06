#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include "pHash.h"
using namespace cimg_library;

void sort_names(char **names, int L1) {
    for (int i = 0; i < L1; i++) {
        int min = i;
        for (int j = i + 1; j < L1; j++) {
            if (strcmp(names[j], names[min]) <= 0) min = j;
        }
        if (i != min) {
            char *swap = names[i];
            names[i] = names[min];
            names[min] = swap;
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 2) {
        printf(" not enough input args\n");
        exit(1);
    }

    const char *dirname1 = argv[1];

    int alpha = 2;
    int level = 1;

    int nbfiles1;
    char **files1 = ph_readfilenames(dirname1, nbfiles1);
    time_t start = time(NULL);
    DP **hashes = ph_mh_imagehashes(files1, nbfiles1);
    printf("Hash computation time: %d\n", time(NULL) - start);
    for (int i = 0; i < nbfiles1; ++i) {
        printf("file: %s\n", hashes[i]->id);
        free(hashes[i]->id);
        free(hashes[i]->hash);
        free(files1[i]);
        free(hashes[i]);
    }
    free(files1);
    free(hashes);
    return 0;
}
