/*

    pHash, the open source perceptual hash library
    Copyright (C) 2009 Aetilius, Inc.
    All rights reserved.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Evan Klinger - eklinger@phash.org
    David Starkweather - dstarkweather@phash.org

*/

#include <getopt.h>
#include <cstdio>
#include "audiophash.h"

int main(int argc, char **argv) {
    struct option long_option[] = {
        {"help", 0, NULL, 'h'},       {"file1", 1, 0, 'f'},
        {"file2", 1, 0, 'g'},         {"threshold", 1, 0, 't'},
        {"block_size", 1, NULL, 'b'}, {NULL, 0, NULL, 0},
    };
    char *file1, *file2;
    float threshold = 0.30;
    int block_size = 256;
    int morehelp = 0;
    while (1) {
        int c;
        if ((c = getopt_long(argc, argv, "f:g:hb:t:", long_option, NULL)) < 0)
            break;
        switch (c) {
            case 'h':
                morehelp++;
                break;
            case 'f':
                file1 = optarg;
                break;
            case 'g':
                file2 = optarg;
                break;
            case 't':
                threshold = atof(optarg);
                break;
            case 'b':
                block_size = atoi(optarg);
        }
    }

    printf("file1: %s\n", file1);
    printf("file2: %s\n", file2);
    printf("threshold: %f\n", threshold);
    printf("blocksize: %d\n", block_size);

    int sr = 8000;
    int channels = 1;

    int N = 0;
    float *buf = ph_readaudio(file1, sr, channels, NULL, N);

    if (!buf) {
        fprintf(stderr, " cannot read file %s, no such file\n", file1);
        return -1;
    }
    int N2 = 0;
    float *buf2 = ph_readaudio(file2, sr, channels, NULL, N2);

    if (!buf2) {
        fprintf(stderr, " cannot read file %s, no such file", file2);
        return -1;
    }

    int nb_frames = 0;
    uint32_t *hash = ph_audiohash(buf, N, sr, nb_frames);

    if (!hash) {
        fprintf(stderr, "unable to calculate hash for %s\n", file1);
        return -1;
    }

    int nb_frames2 = 0;
    uint32_t *hash2 = ph_audiohash(buf2, N2, sr, nb_frames2);

    if (!hash2) {
        fprintf(stderr, "unable to calculate hash for %s\n", file2);
        return -1;
    }

    int Nc = 0;
    double maxC = 0;
    double *C = ph_audio_distance_ber(hash, nb_frames, hash2, nb_frames2,
                                      threshold, block_size, Nc);

    printf("Nc = %d, ", Nc);

    for (int i = 0; i < Nc; i++) {
        if (C[i] > maxC) maxC = C[i];
    }
    printf("cs = %f\n", maxC);

    free(buf);
    free(buf2);
    free(hash);
    free(hash2);
    free(C);

    return 0;
}
