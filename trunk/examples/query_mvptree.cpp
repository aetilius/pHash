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
#include <math.h>
#include "pHash.h"

static int nb_calcs;


float distancefunc(DP *pa, DP *pb){
    nb_calcs++;
    uint8_t *hashA = (uint8_t*)pa->hash;
    uint8_t *hashB = (uint8_t*)pb->hash;
    float d = 10*ph_hammingdistance2(hashA, pa->hash_length, hashB, pb->hash_length);
    float result = exp(d)-1;
    return result;

}


int main(int argc, char **argv){
    if (argc < 3){
	printf("not enough input args\n");
	return 1;
    }

    const char *dir_name = argv[1];/* name of files in directory of query images */
    const char *filename = argv[2];/* name of file to save db */

    float alpha = 2.0f;
    float lvl = 1.0f;

    MVPFile mvpfile;
    ph_mvp_init(&mvpfile);
    mvpfile.filename = strdup(filename);
    mvpfile.hashdist = distancefunc;
    mvpfile.hash_type = BYTEARRAY;
    mvpfile.pgsize = 8192;
    mvpfile.leafcapacity = 46;

    int nbfiles = 0;
    printf("using db %s\n", filename);
    printf("using dir %s for query files\n", dir_name);
    char **files = ph_readfilenames(dir_name,nbfiles);
    if (!files){
	printf("mem alloc error\n");
	exit(1);
    }

    printf("nb query files = %d\n", nbfiles);

    DP *query = ph_malloc_datapoint(mvpfile.hash_type,mvpfile.pathlength);
    
    float radius = 80.0f;
    float threshold = 55.0f;
    int knearest = 20;
    if (argc >= 4){
	char *radius_str = argv[3];
	radius = atof(radius_str);
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
    DP **results = (DP**)malloc(knearest*sizeof(DP**));
    int nbfound = 0, count = 0, sum_calcs = 0;
    int hashlength;
    float d;
    for (int i=0;i<nbfiles;i++){
	printf("query[%d]: %s\n", i, files[i]);

	query->id = files[i];
	query->hash = ph_mh_imagehash(files[i],hashlength,alpha,lvl);
	query->hash_length = hashlength;

	printf("do query ...\n");
	nb_calcs = 0;
	nbfound = 0;
	int res = ph_query_mvptree(&mvpfile,query,knearest,radius,threshold,results,&nbfound);
	if (res != 0){
	    printf("could not complete query, %d\n",res);
	    continue;
	}
        count++;
	sum_calcs += nb_calcs;
	
	printf(" %d files found\n", nbfound);
	for (int i=0;i<nbfound;i++){
	    d = distancefunc(query, results[i]);
	    printf(" %d  %s distance = %f\n", i, results[i]->id, d);
	}
	printf("nb distance calcs: %d\n", nb_calcs);
	free(query->hash);
    } 
   float ave_calcs = (float)sum_calcs/(float)count;      
   printf("ave calcs/query: %f\n", ave_calcs);
    

    return 0;
}
