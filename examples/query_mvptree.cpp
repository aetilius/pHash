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
#include <dirent.h>
#include "pHash.h"

static int nb_calcs;

float imagedistance(DP *pntA, DP *pntB){
    nb_calcs++;
    uint8_t htypeA = pntA->hash_type;
    uint8_t htypeB = pntB->hash_type;
    if (htypeA != htypeB)
	return -1.0;
    if (htypeA != UINT64ARRAY)
	return -1.0;
    if ((pntA->hash_length > 1) || (pntB->hash_length > 1))
	return -1.0;
    ulong64 *hashA = (ulong64*)pntA->hash;
    ulong64 *hashB = (ulong64*)pntB->hash;
    int res = ph_hamming_distance(*hashA, *hashB);
    return (float) res;
}

int main(int argc, char **argv){
 
    const char *dir_name = argv[1];/* name of files in directory of query images */
    const char *filename = argv[2];/* name of file to save db */

    MVPFile mvpfile;
    ph_mvp_init(&mvpfile);
    mvpfile.filename = strdup(filename);
    mvpfile.hashdist = imagedistance;
    mvpfile.hash_type = UINT64ARRAY;

    ulong64 tmphash = 0;

   
    int nbfiles = 0;
    printf("using db %s\n", filename);
    printf("using dir %s for query files\n", dir_name);
    char **files = ph_readfilenames(dir_name,nbfiles);
    if (!files){
	printf("mem alloc error\n");
	exit(1);
    }

    printf("nb query files = %d\n", nbfiles);

    DP *query = NULL;
    float radius = 3.0;
    const int knearest = 20;
    DP **results = (DP**)malloc(knearest * sizeof(DP**));
    int nbfound = 0, count = 0, sum_calcs = 0;

    for (int i=0;i<nbfiles;i++){
	printf("query[%d]: %s ", i, files[i]);
	if (ph_dct_imagehash(files[i],tmphash) < 0){
	    printf("could not get for file %s\n",files[i]);
	    continue;
	}
	printf(" hash = %llx\n", tmphash);
        query = ph_malloc_datapoint(mvpfile.hash_type,mvpfile.pathlength);
        if (!query){
	    printf("mem alloc error\n");
	    continue;
	}
	query->id = strdup(files[i]);
        void *ptr_hash = malloc(8);
	if (!ptr_hash){
            printf("unable to allocate mem\n");
            continue;
        }
        query->hash = ptr_hash;
	ulong64 *ptr = (ulong64*)query->hash;
        *ptr = tmphash;
	query->hash_length = 1;
        count++;
	nb_calcs = 0;
	nbfound = 0;
	int res = ph_query_mvptree(&mvpfile,query,knearest,radius,results,&nbfound);
	if (res != 0){
	    printf("could not complete query\n");
	    continue;
	}
        count++;
	sum_calcs += nb_calcs;

	printf(" %d files found\n", nbfound);
	for (int i=0;i<nbfound;i++){
	    ulong64 *ptrhash = (ulong64*)(results[i]->hash);
	    printf(" %d %s %llx\n", i, results[i]->id, *ptrhash);
	}
	printf("nb distance calcs: %d\n", nb_calcs);
	free(query);
    } 
   float ave_calcs = (float)sum_calcs/(float)count;      
   printf("ave calcs/query: %f\n", ave_calcs);
    

    return 0;
}
