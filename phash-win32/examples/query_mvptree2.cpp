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

#include <math.h>
#include <stdio.h>
#include "pHash.h"

static int nb_calcs;

float distancefunc(DP *pa, DP *pb) {
    nb_calcs++;
    float d = 10 * hammingdistance(pa, pb) / 64;
    float res = exp(d) - 1;
    return res;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        printf("not enough input args\n");
    }
    const char *dir_name =
        argv[1]; /* name of files in directory of query images */
    const char *filename = argv[2]; /* name of file to save db */

    MVPFile mvpfile;
    ph_mvp_init(&mvpfile);
    mvpfile.filename = strdup(filename);
    mvpfile.hashdist = distancefunc;
    mvpfile.hash_type = UINT64ARRAY;
    mvpfile.pgsize = 4096;
    mvpfile.leafcapacity = 25;

    int nbfiles = 0;
    printf("using db %s\n", filename);
    printf("using dir %s for query files\n", dir_name);
    char **files = ph_readfilenames(dir_name, nbfiles);
    if (!files) {
        printf("mem alloc error\n");
        exit(1);
    }

    printf("nb query files = %d\n", nbfiles);

    DP *query = ph_malloc_datapoint(mvpfile.hash_type, mvpfile.pathlength);
    query->hash = malloc(sizeof(ulong64));

    float radius = 50.0f;
    if (argc >= 4) radius = atof(argv[3]);
    int knearest = 20;
    if (argc >= 5) knearest = atoi(argv[4]);
    float threshold = 26.0f;
    if (argc >= 6) threshold = atof(argv[5]);

    DP **results = (DP **)malloc(knearest * sizeof(DP **));
    int nbfound = 0, count = 0, sum_calcs = 0;
    ulong64 tmphash;
    for (int i = 0; i < nbfiles; i++) {
        printf("query[%d]: %s\n", i, files[i]);
        query->id = files[i];
        if (ph_dct_imagehash(files[i], tmphash) < 0) {
            printf("unable to get hash\n");
            continue;
        }
        *((ulong64 *)query->hash) = tmphash;
        query->hash_length = 1;

        printf("do query ...\n");
        nb_calcs = 0;
        nbfound = 0;
        int res = ph_query_mvptree(&mvpfile, query, knearest, radius, threshold,
                                   results, &nbfound);
        if (res != PH_SUCCESS) {
            printf("could not complete query, error %d\n", res);
            continue;
        }
        count++;
        sum_calcs += nb_calcs;

        printf(" %d files found, %d calcs\n", nbfound, nb_calcs);
        for (int i = 0; i < nbfound; i++) {
            printf(" %d  %s distance = %f\n", i, results[i]->id,
                   distancefunc(results[i], query));
        }
        printf("**************************\n");
    }
    float ave_calcs = (float)sum_calcs / (float)count;
    printf("ave calcs/query: %f\n", ave_calcs);

    return 0;
}
