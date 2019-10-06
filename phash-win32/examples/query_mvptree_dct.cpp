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
        printf(
            "usage: progname directory dbname [search radius] [cap] "
            "[threshold]\n");
        return -1;
    }
    const char *dir_name =
        argv[1]; /* name of files in directory of query images */
    const char *filename = argv[2]; /* name of file to save db */

    MVPFile mvpfile;
    mvpfile.filename = strdup(filename);
    mvpfile.hashdist = distancefunc;
    mvpfile.hash_type = UINT64ARRAY;

    int nbfiles = 0;
    printf("using db %s\n", filename);
    printf("using dir %s for query files\n", dir_name);
    char **files = ph_readfilenames(dir_name, nbfiles);
    if (!files) {
        printf("mem alloc error\n");
        return -2;
    }

    printf("nb query files = %d\n", nbfiles);

    DP *query = ph_malloc_datapoint(mvpfile.hash_type);
    if (query == NULL) {
        printf("mem alloc error\n");
        return -3;
    }
    ulong64 *ptmphash = (ulong64 *)malloc(sizeof(ulong64));
    if (ptmphash == NULL) {
        printf("mem alloc error\n");
        return -4;
    }

    float radius = 200.0f;
    if (argc >= 4) radius = atof(argv[3]);
    int knearest = 20;
    if (argc >= 5) knearest = atoi(argv[4]);
    float threshold = 24.0f;
    if (argc >= 6) threshold = atof(argv[5]);

    printf("search radius %f, knearest %d, threshold %f\n", radius, knearest,
           threshold);

    DP **results = (DP **)malloc(knearest * sizeof(DP **));
    if (results == NULL) {
        printf("mem alloc error\n");
        return -5;
    }

    int nbfound = 0, count = 0, sum_calcs = 0;
    for (int i = 0; i < nbfiles; i++) {
        printf("query[%d]: %s\n", i, files[i]);

        if (ph_dct_imagehash(files[i], *ptmphash) < 0) {
            printf("unable to get hash\n");
            continue;
        }

        query->id = files[i];
        query->hash = (void *)ptmphash;
        query->hash_length = 1;

        nb_calcs = 0;
        nbfound = 0;
        int res = ph_query_mvptree(&mvpfile, query, knearest, radius, threshold,
                                   results, &nbfound);
        if (res != PH_SUCCESS && res != PH_ERRCAP) {
            printf("could not complete query, error %d\n", res);
            continue;
        }
        count++;
        sum_calcs += nb_calcs;

        printf(" %d files found in %d distance calculations\n", nbfound,
               nb_calcs);

        nbfound = (nbfound < knearest) ? nbfound : knearest;
        for (int j = 0; j < nbfound; j++) {
            printf(" =====> result[%d] %s distance = %f\n", j, results[j]->id,
                   distancefunc(results[j], query));
            hfree(results[j]->id);
            hfree(results[j]->hash);
            ph_free_datapoint(results[j]);
        }

        printf("\n\n");
    }
    float ave_calcs = (float)sum_calcs / (float)count;
    printf("\n\n\n\nave calcs/query: %f\n", ave_calcs);

    return 0;
}
