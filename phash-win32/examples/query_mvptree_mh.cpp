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

static int nb_calcs;

float distancefunc(DP *pa, DP *pb) {
    nb_calcs++;
    uint8_t *hashA = (uint8_t *)pa->hash;
    uint8_t *hashB = (uint8_t *)pb->hash;
    float d = 10 * ph_hammingdistance2(hashA, pa->hash_length, hashB,
                                       pb->hash_length);
    float res = exp(d) - 1;
    return res;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        printf("not enough input args\n");
        return -1;
    }

    const char *dir_name =
        argv[1]; /* name of files in directory of query images */
    const char *filename = argv[2]; /* name of file to save db */

    float alpha = 2.0f;
    float lvl = 1.0f;

    MVPFile mvpfile;
    mvpfile.filename = strdup(filename);
    mvpfile.hashdist = distancefunc;
    mvpfile.hash_type = BYTEARRAY;

    int nbfiles = 0;
    printf("using db %s\n", filename);
    printf("using dir %s for query files\n", dir_name);
    char **files = ph_readfilenames(dir_name, nbfiles);
    if (!files) {
        printf("mem alloc error\n");
        return -2;
    }

    printf("nb query files = %d\n", nbfiles);

    float radius = 200.0f;
    if (argc >= 4) {
        radius = atof(argv[3]);
    }
    int knearest = 2;
    if (argc >= 5) {
        knearest = atoi(argv[4]);
    }
    float threshold = 55.0f;
    if (argc >= 6) {
        threshold = atof(argv[5]);
    }
    printf("radius = %f, knearest = %d, threshold = %f\n", radius, knearest,
           threshold);

    DP *query = ph_malloc_datapoint(mvpfile.hash_type);

    DP **results = (DP **)malloc(knearest * sizeof(DP **));
    if (results == NULL) {
        printf("unable to allocate memory\n");
        return -3;
    }

    int nbfound = 0, count = 0, sum_calcs = 0;
    int hashlength;
    for (int i = 0; i < nbfiles; i++) {
        printf("query[%d]: %s\n", i, files[i]);
        uint8_t *hasht = ph_mh_imagehash(files[i], hashlength, alpha, lvl);
        if (hasht == NULL) {
            printf("unable to get hash\n");
            continue;
        }

        query->id = files[i];
        query->hash = hasht;
        query->hash_length = hashlength;

        nb_calcs = 0;
        nbfound = 0;
        int res = ph_query_mvptree(&mvpfile, query, knearest, radius, threshold,
                                   results, &nbfound);
        if (res != PH_SUCCESS && res != PH_ERRCAP) {
            printf("unable to complete query - %d\n", res);
            hfree(query->hash);
            continue;
        }
        printf("\n\n%d found in %d calculations\n", nbfound, nb_calcs);

        count++;
        sum_calcs += nb_calcs;

        for (int j = 0; j < nbfound; j++) {
            float d = distancefunc(query, results[j]);
            printf("====>file: %s dist = %f\n", results[j]->id, d);
            hfree(results[j]->id);
            hfree(results[j]->hash);
            hfree(results[j]);
        }
        printf("\n\n");
        hfree(query->hash);
    }
    float ave_calcs = (float)sum_calcs / (float)count;
    printf("\n\nave calcs/query: %f\n", ave_calcs);

    return 0;
}
