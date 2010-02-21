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

    const int sr = 8000;
    const int nbchannels = 1;

    MVPFile mvpfile;
    ph_mvp_init(&mvpfile);
    mvpfile.filename = strdup(filename);
    mvpfile.hashdist = audiohashdistance;
    mvpfile.hash_type = UINT32ARRAY;

    int nbfiles = 0;
    printf("dir name: %s\n", dir_name);
    char **files = ph_readfilenames(dir_name,nbfiles);
    if (!files){
	printf("mem alloc error\n");
	exit(1);
    }
    printf("nbfiles = %d\n", nbfiles);
    DP **hashlist = (DP**)malloc(nbfiles*sizeof(DP*));
    if (!hashlist){
	printf("mem alloc error\n");
	exit(1);
    }

    int count = 0;
    float *buf;
    int N = 0;
    uint32_t *hash;
    int hashlen = 0;
    DP *dp;
    for (int i=0;i<nbfiles;i++){
	printf("file[%d]: %s ", i, files[i]);
	buf = ph_readaudio(files[i], sr, nbchannels, N);
	if (!buf){
	    printf("unable to read file\n");
	    continue;
	}
	printf("N = %d ", N);

	hash = ph_audiohash(buf, N, sr, hashlen);
	if (!hash){
	    printf("unable to get hash\n");
	    free(buf);
	}
	printf("nb hashes = %d\n", hashlen);

        dp = ph_malloc_datapoint(mvpfile.hash_type,mvpfile.pathlength);
	if (!dp){
	    printf("mem alloc error\n");
	    free(buf);
	    free(hash);
	}
	dp->id = files[i];
	dp->hash = (void*)hash;
	dp->hash_length = hashlen;

	hashlist[count++] = dp;
    }

    printf("add %d files to file %s\n", count, filename);
    int nbsaved;
    MVPRetCode ret = ph_add_mvptree(&mvpfile, hashlist, count,nbsaved);
    if (ret != PH_SUCCESS){
		printf("unable to add points to %s\n", filename);
    }
    printf("save %d files out of %d\n", nbsaved, count);

    free(hashlist);

    return 0;
}
