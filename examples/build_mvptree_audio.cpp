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

float audiohashdistance(DP *dpA, DP *dpB){

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
    return (float)res;

}

int main(int argc, char **argv){
 
    const char *dir_name = argv[1];/* name of dir to retrieve image files */
    const char *filename = argv[2];/* name of file to save db */

    MVPFile mvpfile;
    mvpfile.branchfactor = 2;
    mvpfile.leafcapacity = 25;
    mvpfile.pathlength = 5;
    mvpfile.pgsize = (1 << 18);
    mvpfile.filename = strdup(filename);
    mvpfile.hashdist = audiohashdistance;
    mvpfile.hash_type = UINT32ARRAY;

    int sr = 8000;
    int nbchannels = 1;

    int nbfiles = 0;
    printf("dir name: %s\n", dir_name);
    char **files = ph_readfilenames(dir_name,nbfiles);
    if (!files){
	printf("mem alloc error\n");
	exit(1);
    }
    printf("nbfiles = %d\n", nbfiles);
    DP **hashlist = (DP**)malloc(nbfiles*sizeof(DP**));
    if (!hashlist){
	printf("mem alloc error\n");
	exit(1);
    }

    int count = 0;
    float *buf = NULL;
    uint32_t *hash = NULL;
    DP *dp = NULL;
    int nbsamples = 0;
    int hashlength = 0;
    for (int i=0;i<nbfiles;i++){
	buf = ph_readaudio(files[i], sr, nbchannels, nbsamples);
	if (!buf){
	    printf("unable to read file: %s\n", files[i]);
	    continue;
	}
	hash = ph_audiohash(buf, nbsamples, sr, hashlength);
	if (!hash){
	    printf("unable to compute hash: %s\n", files[i]);
	    continue;
	}
	printf("file[%d]: %s, samples %d, hashlen = %d\n", i, files[i],nbsamples,hashlength);
        dp = ph_malloc_datapoint(mvpfile.hash_type,mvpfile.pathlength);
        if (!dp){
	    printf("mem alloc error\n");
	    exit(1);
	}
	dp->id = files[i];
	dp->hash = hash;
	dp->hash_length = hashlength;
	hashlist[count++] = dp;

	free(buf);
    }

    printf("save %d files\n", count);

    if (ph_save_mvptree(&mvpfile, hashlist, count) < 0){
	printf("unable to save %s\n", filename);
	exit(1);
    }

    return 0;
}
