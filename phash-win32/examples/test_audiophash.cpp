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

#include <stdio.h>
#include <string.h>
#include <vector>
#include "audiophash.h"
#include "phash.h"

#define TRUE 1
#define FALSE 0

using namespace std;

typedef struct hashpoint {
    char *filename;
    uint32_t *hash;
    int length;
} HP;

int main(int argc, char **argv) {
    const char *msg = ph_about();
    printf(msg);
    printf("\n\n");

    if (argc < 3) {
        printf("no args");
        printf("usage: %s dir1 dir2\n", argv[0]);
        return -1;
    }
    const float nbsecs = 45.0f;
    const char *dir_name = argv[1];   // first directory
    const char *dir_name2 = argv[2];  // second directory
    const float threshold = 0.30;     // ber threshold (0.25-0.35)
    const int block_size = 256;       // number of frames to compare at a time
    const int sr = 8000;              // sample rate to convert the stream
    const int channels = 1;           // number of channels to convert stream

    int N1;
    char **files1 = ph_readfilenames(dir_name, N1);
    int N2;
    char **files2 = ph_readfilenames(dir_name2, N2);
    printf("dir1: %s %d\n", dir_name, N1);
    printf("dir2: %s %d\n", dir_name2, N2);

    if (N1 != N2) {
        printf("unequal number files in both directories\n");
        return -1;
    }
    vector<HP> HashList(N1);

    float *tmpbuf;
    float *buf = new float[1 << 25];
    int buflen = (1 << 25);
    if (!buf) {
        printf("unable to alloc buffer storage\n");
        return -1;
    }

    int hashN = (1 << 26);
    uint32_t *hashes = new uint32_t[1 << 26];
    if (!hashes) {
        printf("unable to alloc storage for hashes\n");
        return -1;
    }

    int hash_index = 0;
    uint32_t *hash1, *hash2;
    int buflen1, buflen2;
    int nbframes1, nbframes2;
    double *ptrC;
    int Nc;
    printf("intra distances\n");
    for (int i = 0; i < N1; i++) {
        printf("  files1[%d] = %s\n", i, files1[i]);
        buflen1 = buflen;
        tmpbuf = ph_readaudio(files1[i], sr, channels, buf, buflen1, nbsecs);
        if (!tmpbuf) continue;
        buf = tmpbuf;

        hash1 = ph_audiohash(buf, buflen1, hashes, hashN, sr, nbframes1);
        if (!hash1) continue;
        hashes += nbframes1;
        hashN -= nbframes1;

        HashList[i].filename = files1[i];
        HashList[i].hash = hash1;
        HashList[i].length = nbframes1;

        printf("  files2[%d] = %s\n", i, files2[i]);
        buflen2 = buflen;
        tmpbuf = ph_readaudio(files2[i], sr, channels, buf, buflen2, nbsecs);
        if (!tmpbuf) continue;
        buf = tmpbuf;

        hash2 = ph_audiohash(buf, buflen2, hashes, hashN, sr, nbframes2);
        if (!hash2) continue;
        hashes += nbframes2;
        hashN -= nbframes2;

        ptrC = ph_audio_distance_ber(hash1, nbframes1, hash2, nbframes2,
                                     threshold, block_size, Nc);

        double maxC = 0.0f;
        for (int j = 0; j < Nc; j++) {
            if (ptrC[j] > maxC) maxC = ptrC[j];
        }
        printf("distance %f\n", maxC);
        printf("********************************\n");
        delete[] ptrC;
    }

    printf("hit any key to continue\n");
    getchar();
    for (int i = 0; i < N1; i++) {
        for (int j = i + 1; j < N1; j++) {
            ptrC = ph_audio_distance_ber(HashList[i].hash, HashList[i].length,
                                         HashList[j].hash, HashList[j].length,
                                         threshold, block_size, Nc);
            double maxC = 0.0f;
            for (int k = 0; k < Nc; k++) {
                if (ptrC[k] > maxC) maxC = ptrC[k];
            }
            printf(" %d %d dist = %f\n", i, j, maxC);
            free(ptrC);
        }
    }

    return 0;
}
