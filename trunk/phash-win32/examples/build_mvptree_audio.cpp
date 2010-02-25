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
#include "audiophash.h"


float distancefunc(DP *pa, DP *pb){
  int Nc;
  float threshold =  0.30f;
  int block_size = 256;
  double *ptrC = ph_audio_distance_ber((uint32_t*)pa->hash,(int)pa->hash_length,(uint32_t*)pb->hash,(int)pb->hash_length,threshold,block_size,Nc);

  double maxC = 0.0;
  for (int i=0;i<Nc;i++){
	  if (ptrC[i] > maxC){
          maxC = ptrC[i];
	  }
  }  
  if (ptrC) 
	  free(ptrC);

  double d = 10*(1.0 - maxC);
  float res = (float)(exp(d) - 1.0);
  return res;

}

int main(int argc, char **argv){
	if (argc < 3){
       printf("not enough input args\n");
       return 1;
	}
    const char *dir_name = argv[1];/* name of dir to retrieve image files */
    const char *filename = argv[2];/* name of indexing db, e.g. 'audiodb'  */

    const int sr = 8000;            /* convert to sr */
    const int nbchannels = 1;       /* convert to number of channels */
    const float nbsecs = 45.0f;     /* nb secs of audio to read from each file */

    MVPFile mvpfile;                /* mvp tree indexing configuration */ 
    mvpfile.branchfactor = 2;
    mvpfile.pathlength = 5;
    mvpfile.leafcapacity = 44;
    mvpfile.pgsize = 1 << 19; /* 2^29=524,288 */
    mvpfile.filename = strdup(filename);
    mvpfile.hashdist = distancefunc;
    mvpfile.hash_type = UINT32ARRAY;

    int nbfiles = 0;
    printf("dir name: %s\n", dir_name);
    char **files = ph_readfilenames(dir_name,nbfiles);
    if (!files){
		printf("mem alloc error\n");
		return 1;
    }
    printf("nbfiles = %d\n", nbfiles);
    DP **hashlist = (DP**)malloc(nbfiles*sizeof(DP*));
    if (!hashlist){
	    printf("mem alloc error\n");
	    return 1;
    }
    
    /* audio buffer to pass to readaudio() */ 
    int count = 0;
    float *sigbuf = (float*)malloc(1<<21);
    int buflen = (1<<21)/sizeof(float);
    int N;
  
    /*buffer for all hashes */ 
    uint32_t *hashes = (uint32_t*)malloc(1<<21);
    int hashspacelength = (1<<21)/sizeof(uint32_t);
    int hashspaceleft = hashspacelength;
  
    uint32_t *hash = hashes;
    int nbframes;
    for (int i=0;i<nbfiles;i++){
		printf("file[%d]: %s \n", i, files[i]);
        hashlist[count] = ph_malloc_datapoint(mvpfile.hash_type);
		if (hashlist[count] == NULL){
			printf("mem alloc error\n");
			continue;
		}
		hashlist[count]->id = files[i];
        N = buflen;
        float *buf = ph_readaudio(files[i], sr, nbchannels, sigbuf, N, nbsecs);
		if (buf == NULL){
            printf("unable to get signal\n");
            ph_free_datapoint(hashlist[count]);
            continue;
		}
        printf("nb sampels = %d\n", N);
        uint32_t *hash1 = ph_audiohash(buf, N, hash, hashspaceleft, sr, nbframes);
		if (hash1 == NULL){
			printf("unable to get hash\n\n");
            ph_free_datapoint(hashlist[count]);
		    continue;
		}
        printf("nb hashes %d\n", nbframes);
        hash += nbframes;
        hashspaceleft -= nbframes;

        hashlist[count]->hash = hash1;
		hashlist[count]->hash_length = nbframes;
        count++;
    }
    printf("save files into tree ...\n");
    MVPRetCode errcode = ph_save_mvptree(&mvpfile, hashlist, count);
    if (errcode != PH_SUCCESS){
		printf("unable to save %s, err %d\n", filename, errcode);
		return 1;
    }
    printf("saved files\n");
	
    for (int i=0;i<nbfiles;i++){
        free(files[i]);
        free(hashlist[i]);
	}
    free(files);
    free(hashlist);
    free(sigbuf);
    return 0;
}
