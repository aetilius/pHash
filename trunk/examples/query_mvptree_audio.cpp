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
#include "audiophash.h"

static int nb_calcs;

float audiohashdistance(DP *dpA, DP *dpB){

    nb_calcs++;

    uint32_t *hash1 = (uint32_t*)dpA->hash;
    int N1 = dpA->hash_length;
    uint32_t *hash2 = (uint32_t*)dpB->hash;
    int N2 = dpB->hash_length;

    float threshold = 0.30;
    int blocksize = 256;
    int Nc=0;
    double *ptrC = ph_audio_distance_ber(hash1, N1, hash2, N2, threshold,blocksize,Nc);

    double maxC = 0;
    for (int i=0;i<Nc;i++){
	if (ptrC[i] > maxC)
	    maxC = ptrC[i];
    }
    double res = 1000*(1-maxC);
    free(ptrC);
    return (float)res;
}

int main(int argc, char **argv){
    if (argc < 3){
	printf("not enough input args\n");
        printf("usage: %s directory dbname radius knearest threshold\n", argv[0]);
        return -1;
    } 
 
    const char *dir_name = argv[1];/* name of files in directory of query images */
    const char *filename = argv[2];/* name of file to save db */

    const int sr = 8000;
    const int nbchannels = 1;
    const float nbsecs = 45.0f;

    MVPFile mvpfile;
    mvpfile.filename = strdup(filename);
    mvpfile.hashdist = audiohashdistance;
    mvpfile.hash_type = UINT32ARRAY;
   
    int nbfiles = 0;
    printf("using db %s\n", filename);
    printf("using dir %s for query files\n", dir_name);
    char **files = ph_readfilenames(dir_name,nbfiles);
    if (!files){
	printf("mem alloc error\n");
	return -2;
    }
    printf("nb files = %d\n", nbfiles);

    float radius = 700.0;
    if (argc >= 3){
	radius = atof(argv[3]);
    }
    int knearest = 20;
    if (argc >= 4){
	knearest = atoi(argv[4]);
    }
    float threshold = 300.0;
    if (argc >= 5){
        threshold = atof(argv[5]);
    }
    printf("radius %f, knearest %d, threshold %f\n", radius, knearest, threshold);

    DP **results = (DP**)malloc(knearest*sizeof(DP*));
    if (!results){
	printf("mem alloc error\n");
	return -3;
    }
    DP *query = ph_malloc_datapoint(mvpfile.hash_type);

    int nbfound = 0, count = 0, sum_calcs = 0;
    float *sigbuf = (float*)malloc(1<<21);
    int buflen = (1<<21)/sizeof(float);
    for (int i=0;i<nbfiles;i++){
        printf("query[%d]: %s\n", i, files[i]);
        int N = buflen;
	float *buf = ph_readaudio(files[i], sr, nbchannels, sigbuf, N, nbsecs);
	if (!buf){
	    printf("could not read audio\n");
	    continue;
	}
        printf("samples %d\n", N);
        int nbframes = 0;
	uint32_t *hash = ph_audiohash(buf,N,sr,nbframes);
	if (!hash){
	    printf("could not get hash\n");
	    continue;
	}
	printf("number hash frames %d\n\n", nbframes);

	query->id = files[i];
	query->hash = hash;
	query->hash_length = nbframes;

	count++;
	nb_calcs = 0;
	nbfound = 0;
	int ret = ph_query_mvptree(&mvpfile,query,knearest,radius,threshold, results,nbfound);
	if (ret != PH_SUCCESS && ret != PH_ERRCAP){
	    printf("could not complete query - retcode %d\n",ret);
            free(hash);
	    continue;
	}
	sum_calcs += nb_calcs;

	printf("search: %d found,  %d distance calc's\n", nbfound,nb_calcs);
	for (int j=0;j<nbfound;j++){
	    float dist = audiohashdistance(query,results[j]);
	    printf("     %d  %s distance = %f\n", j, results[j]->id,dist);
	}
        free(hash);
        if (sigbuf != buf){
	    free(buf);
	}
        for (int j=0;j<nbfound;j++){
	    free(results[j]->id);
            free(results[j]->hash);
            ph_free_datapoint(results[j]);
	}
	printf("\n\n**********************************************\n");
    } 
   float ave_calcs = (float)sum_calcs/(float)count;      
   printf("ave calcs/query: %f\n", ave_calcs);
    
    return 0;
}
