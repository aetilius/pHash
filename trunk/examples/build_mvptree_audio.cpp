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
#include "audiophash.h"

float audiohashdistance(DP *dpA, DP *dpB){

    uint32_t *hash1 = (uint32_t*)dpA->hash;
    int N1 = dpA->hash_length;
    uint32_t *hash2 = (uint32_t*)dpB->hash;
    int N2 = dpB->hash_length;

    float threshold = 0.30;
    int blocksize = 256;
    int Nc=0;
    double *ptrC = ph_audio_distance_ber(hash1, N1, hash2, N2, threshold,blocksize,Nc);

    double maxC = 0.0f;
    for (int i=0;i<Nc;i++){
	if (ptrC[i] > maxC)
	    maxC = ptrC[i];
    }
    free(ptrC);
    double res = 1000*(1-maxC);
    return (float)res;

}

int main(int argc, char **argv){
 
    const char *dir_name = argv[1];/* name of dir to retrieve image files */
    const char *filename = argv[2];/* name of file to save db */

    MVPFile mvpfile;
    mvpfile.branchfactor = 2;
    mvpfile.leafcapacity = 40;
    mvpfile.pathlength = 5;
    mvpfile.pgsize = (1 << 19); /* 524,288 */ 
    mvpfile.filename = strdup(filename);
    mvpfile.hashdist = audiohashdistance;
    mvpfile.hash_type = UINT32ARRAY;

    int sr = 8000;
    int nbchannels = 1;
    float nbsecs = 45.0f;

    int nbfiles = 0;
    printf("dir name: %s\n", dir_name);
    char **files = ph_readfilenames(dir_name,nbfiles);
    if (!files){
	printf("unable to read files from directory %s\n",dir_name);
	return -2;
    }
    printf("nbfiles = %d\n", nbfiles);
    DP **hashlist = (DP**)malloc(nbfiles*sizeof(DP**));
    if (!hashlist){
	printf("mem alloc error\n");
	return -3;
    }

    int count = 0;
    float *sigbuf = (float*)malloc(1<<21);
    int buflen = (1<<21)/sizeof(float);

    for (int i=0;i<nbfiles;i++){
        printf("file[%d]: %s\n", i, files[i]);
        int N = buflen;
	float *buf = ph_readaudio(files[i], sr, nbchannels, sigbuf, N, nbsecs);
	if (!buf){
	    printf("unable to read file: %s\n", files[i]);
	    continue;
	}
	printf("signal length %d\n", N);
        int nbframes=0;
	uint32_t *hash = ph_audiohash(buf, N, sr, nbframes);
	if (!hash){
	    printf("unable to compute hash\n");
	    continue;
	}
	printf("nbframes %d\n", nbframes);
        hashlist[count] = ph_malloc_datapoint(mvpfile.hash_type);
        if (!hashlist[count]){
	    printf("mem alloc error\n");
            free(hash);
	    return -4;
	}

	hashlist[count]->id = files[i];
        hashlist[count]->hash = hash;
        hashlist[count++]->hash_length = nbframes;

	if (buf != sigbuf)
            free(buf);
    }

    printf("save %d files\n", count);

    MVPRetCode ret = ph_save_mvptree(&mvpfile, hashlist, count);
    if (ret != PH_SUCCESS){
	printf("unable to save %s file - ret code %d\n", filename, ret);
    }
    printf("done\n");

    return 0;
}
