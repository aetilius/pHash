#include <stdio.h>
#include <stdlib.h>
#include "pHash.h"

int main(int argc, char **argv) {
    if (argc < 3) {
        printf("not enough input args\n");
        exit(1);
    }
    const char *file1 = argv[1];
    const char *file2 = argv[2];

    printf("file1: %s\n", file1);
    int nbhashes1 = 0;
    TxtHashPoint *hash1 = ph_texthash(file1, &nbhashes1);
    if (!hash1) {
        printf("unable to complete text hash function\n");
        exit(1);
    }
    printf("length %d\n", nbhashes1);

    printf("file2: %s\n", file2);
    int nbhashes2 = 0;
    TxtHashPoint *hash2 = ph_texthash(file2, &nbhashes2);
    if (!hash2) {
        printf("unable to complete text hash function\n");
        exit(2);
    }

    printf("length %d\n", nbhashes2);
    int count, j;
    TxtMatch *matches =
        ph_compare_text_hashes(hash1, nbhashes1, hash2, nbhashes2, &count);
    if (!matches) {
        printf("unable to complete compare function\n");
        exit(3);
    }

    printf(" %d matches\n", count);
    printf(" indxA  indxB  length\n");
    for (j = 0; j < count; j++) {
        printf(" %d %d %d\n", matches[j].first_index, matches[j].second_index,
               matches[j].length);
    }

    return 0;
}
