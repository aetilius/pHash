#include <stdio.h>
#include <stdlib.h>
#include "pHash.h"

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
    if (argc < 3) {
        printf("not enough input args\n");
        exit(1);
    }
    const char *dir1 = argv[1];
    const char *dir2 = argv[2];

    int N1 = 0;
    char **files1 = ph_readfilenames(dir1, N1);
    sort_names(files1, N1);
    int N2 = 0;
    char **files2 = ph_readfilenames(dir2, N2);
    sort_names(files2, N2);

    if (N1 != N2) {
        printf("unequal number files in each directory\n");
        exit(1);
    }

    TxtHashPoint *hash1, *hash2;
    TxtMatch *matches;
    int nbhashes1, nbhashes2, count;
    for (int i = 0; i < N1; i++) {
        hash1 = ph_texthash(files1[i], &nbhashes1);
        if (!hash1) {
            printf("unable to get hash\n");
            continue;
        }
        printf("file%d: %s, length %d\n", i, files1[i], nbhashes1);
        hash2 = ph_texthash(files2[i], &nbhashes2);
        if (!hash2) {
            printf("unable to get hash\n");
            free(hash1);
            continue;
        }
        printf("file%d: %s, length %d\n", i, files2[i], nbhashes2);

        matches =
            ph_compare_text_hashes(hash1, nbhashes1, hash2, nbhashes2, &count);
        if (!matches) {
            printf("no matches found\n");
            free(hash1);
            free(hash2);
            continue;
        }
        int maxlength = 0;
        for (int j = 0; j < count; j++) {
            if (matches[j].length > maxlength) maxlength = matches[j].length;
        }
        printf("     no. matches %d, max length %d\n", count, maxlength);

        free(hash1);
        free(hash2);
        free(matches);
    }
    getchar();
    printf("******************\n");
    for (int i = 0; i < N1; i++) {
        hash1 = ph_texthash(files1[i], &nbhashes1);
        if (!hash1) continue;

        for (int j = i + 1; j < N1; j++) {
            hash2 = ph_texthash(files1[j], &nbhashes2);
            if (!hash2) continue;
            matches = ph_compare_text_hashes(hash1, nbhashes1, hash2, nbhashes2,
                                             &count);
            if (!matches) continue;
            printf("file%d: %s\n", i, files1[i]);
            printf("file%d: %s\n", j, files2[j]);
            int maxlength = 0;
            for (int k = 0; k < count; k++) {
                if (matches[k].length > maxlength)
                    maxlength = matches[k].length;
            }
            printf("  no. matches %d, max length %d\n\n", count, maxlength);

            free(hash2);
            free(matches);
        }
        free(hash1);
    }

    printf("done\n");
    return 0;
}
