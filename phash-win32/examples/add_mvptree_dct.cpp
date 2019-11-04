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
#include "pHash.h"

float distancefunc(DP *pa, DP *pb) {
    float d = 10 * hammingdistance(pa, pb) / 64;
    float res = exp(d) - 1;
    return res;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        printf("not enough input args\n\n");
        printf("usage: progname directory dbname\n");
        return -1;
    }

    const char *dir_name = argv[1]; /* name of dir to retrieve image files */
    const char *filename = argv[2]; /* name of file to save db */

    MVPFile mvpfile;
    mvpfile.filename = strdup(filename);
    mvpfile.hashdist = distancefunc;
    mvpfile.hash_type = UINT64ARRAY;

    int nbfiles = 0;
    printf("dir name: %s\n", dir_name);
    char **files = ph_readfilenames(dir_name, nbfiles);
    if (!files) {
        printf("mem alloc error\n");
        return -2;
    }
    printf("nbfiles = %d\n", nbfiles);
    DP **hashlist = (DP **)malloc(nbfiles * sizeof(DP *));
    if (!hashlist) {
        printf("mem alloc error\n");
        return -3;
    }
    ulong64 tmphash;
    int count = 0;
    for (int i = 0; i < nbfiles; i++) {
        printf("file[%d]: %s\n", i, files[i]);
        hashlist[count] = ph_malloc_datapoint(mvpfile.hash_type);
        if (hashlist[count] == NULL) {
            printf("mem alloc error\n");
            return -4;
        }
        hashlist[count]->hash = malloc(sizeof(ulong64));
        if (hashlist[count]->hash == NULL) {
            printf("mem alloc error\n");
            ph_free_datapoint(hashlist[count]);
            return -5;
        }
        if (ph_dct_imagehash(files[i], tmphash) < 0) {
            printf("unable to get hash\n");
            hfree(hashlist[count]->hash);
            ph_free_datapoint(hashlist[count]);
            continue;
        }
        hashlist[count]->id = files[i];
        *((ulong64 *)hashlist[count]->hash) = tmphash;
        hashlist[count]->hash_length = 1;
        count++;
    }
    int nbsaved;
    printf("add files to file %s\n", filename);
    MVPRetCode ret = ph_add_mvptree(&mvpfile, hashlist, count, nbsaved);
    if (ret != PH_SUCCESS) {
        printf("unable to add files\n");
        return -6;
    }
    printf("added %d of %d files to %s\n", nbsaved, count, mvpfile.filename);

    return 0;
}
