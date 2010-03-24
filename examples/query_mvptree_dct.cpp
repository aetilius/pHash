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
    D Grant Starkweather - dstarkweather@phash.org

*/

#include <stdio.h>
#include <math.h>
#include "pHash.h"

static int nb_calcs;

float distancefunc(DP *pa, DP *pb){
    nb_calcs++;
    float d = ph_hamming_distance(*((ulong64*)pa->hash),*((ulong64*)pb->hash));
    return d;
}

int main(int argc, char **argv){
    if (argc < 3){
	printf("not enough input args\n");
        printf("usage: %s directory dbname [radius] [knearest] [threshold]\n", argv[0]);
	return -1;
    }

    const char *dir_name = argv[1];/* name of files in directory of query images */
    const char *filename = argv[2];/* name of file to save db */

    MVPFile mvpfile;
    mvpfile.filename = strdup(filename);
    mvpfile.hashdist = distancefunc;
    mvpfile.hash_type = UINT64ARRAY;

    int nbfiles = 0;
    printf("using db %s\n", filename);
    printf("using dir %s for query files\n", dir_name);
    char **files = ph_readfilenames(dir_name,nbfiles);
    if (!files){
	printf("unable to read files from directory\n");
	return -2;
    }

    printf("nb query files = %d\n", nbfiles);

    DP *query = ph_malloc_datapoint(mvpfile.hash_type);
    if (query == NULL){
        printf("mem alloc error\n");
        return -3;
    }

    float radius = 30.0f;
    float threshold = 15.0f;
    int knearest = 20;
    if (argc >= 4){
	radius = atof(argv[3]);
    }
    if (argc >= 5){
	knearest = atoi(argv[4]);
    }
    if (argc >= 6){
        threshold = atof(argv[5]);
    }
    printf("radius = %f\n", radius);
    printf("knearest = %d\n", knearest);
    printf("threshold = %f\n", threshold);

    DP **results = (DP**)malloc(knearest*sizeof(DP*));
    if (results == NULL){
        return -3;
    }
    ulong64 tmphash = 0x0000000000000000;
    int nbfound = 0, count = 0, sum_calcs = 0;
    for (int i=0;i<nbfiles;i++){

        if (ph_dct_imagehash(files[i],tmphash) < 0){
	    printf("unable to get hash\n");
            continue;
	}
	printf("query[%d]: %s %llx\n", i, files[i], tmphash);
        query->id = files[i];
        query->hash = &tmphash;
        query->hash_length = 1;

	nb_calcs = 0;
	nbfound = 0;
	MVPRetCode retcode = ph_query_mvptree(&mvpfile,query,knearest,radius,threshold,results,nbfound);
	if (retcode != PH_SUCCESS && retcode != PH_ERRCAP){
	    printf("could not complete query, %d\n",retcode);
	    continue;
	}
        count++;
	sum_calcs += nb_calcs;
	
	printf(" %d files found\n", nbfound);
	for (int j=0;j<nbfound;j++){
	    float d = distancefunc(query, results[j]);
	    printf(" %d  %s distance = %f\n", j, results[j]->id, d);
	}
	printf("nb distance calcs: %d\n", nb_calcs);
	for (int j=0;j<nbfound;j++){
            free(results[j]->id);
            results[j]->id = NULL;
            free(results[j]->hash);
            results[j]->id = NULL;
	    ph_free_datapoint(results[j]);
	}

    }
 
   float ave_calcs = (float)sum_calcs/(float)count;      
   printf("ave calcs/query: %f\n", ave_calcs);
    

   for (int i=0;i<nbfiles;i++){
       free(files[i]);
   }
   free(files);

   ph_free_datapoint(query);
   free(results);
   free(mvpfile.filename);

    return 0;
}

