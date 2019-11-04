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

#include "pHash.h"
#include "stdio.h"

int main(int argc, char **argv) {
    if (argc < 2) {
        printf("Not enough arguments\n");
        printf(
            "Please specify a list of video paths. Paths after the first will "
            "be compared to the first.\n");
        return 1;
    }

    ulong64 **hashes = (ulong64 **)malloc(sizeof(ulong64 *) * (argc - 1));
    int *lengths = (int *)malloc(sizeof(int) * (argc - 1));

    printf("Comparing additional paths to: %s\n\n", argv[1]);

    for (int i = 1; i < argc; i++) {
        printf("Hashing %s...\n", argv[i]);
        hashes[i - 1] = ph_dct_videohash(argv[i], lengths[i - 1]);
        if (hashes[i - 1] == NULL) {
            printf("Failed to hash video: %s\n", argv[i]);
            break;
        }
        printf("%s: %lx (length %d)\n", argv[i], hashes[i - 1][0],
               lengths[i - 1]);
        // If this isn't the first hash, then compare to the first.
        if (i > 1) {
            float similarity = ph_dct_videohash_dist(
                hashes[0], lengths[0], hashes[i - 1], lengths[i - 1]);
            printf("%s similarity: %f\n", argv[i], similarity);
        }
        printf("\n");
    }

    int result = 0;
    for (int i = 0; i < argc - 1; i++) {
        if (hashes[i] == NULL) {
            result = -1;
        } else {
            free(hashes[i]);
        }
    }
    free(hashes);
    free(lengths);

    return result;
}
