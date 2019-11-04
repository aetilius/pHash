#include <stdio.h>
#include "audiophash.h"
#include "pHash.h"

int main(int argc, char **argv) {
    if (argc < 2) {
        printf("not enough input args\n");
        return 1;
    }

    const char *dirname = argv[1];

    const float nbsecs = 0.0;
    const int sr = 8000;
    const int chs = 1;

    float *sigbuf = (float *)malloc(1 << 26);
    const int N = (1 << 26) / sizeof(float);

    int nbfiles;
    char **files = ph_readfilenames(dirname, nbfiles);
    if (files == NULL) {
        printf("unable to read filenames in dir\n");
        return 2;
    }

    printf("nbfiles %d\n", nbfiles);
    for (int i = 0; i < nbfiles; i++) {
        printf("file[%d]: %s\n", i, files[i]);
        int buflen = N;
        float *buf = ph_readaudio(files[i], sr, chs, sigbuf, buflen, nbsecs);
        if (buf == NULL) {
            printf("unable to extract audio\n");
            continue;
        }
        printf("buf = %p, len = %d\n", buf, buflen);

        int nbframes;
        uint32_t *hash = ph_audiohash(buf, buflen, sr, nbframes);

        if (hash == NULL) {
            printf("unable to calculate hash\n");
            continue;
        }

        printf("hash frames %d\n", nbframes);

        free(hash);
        if (buf != sigbuf) free(buf);
    }

    free(sigbuf);
    for (int i = 0; i < nbfiles; i++) {
        free(files[i]);
    }
    free(files);

    return 0;
}
