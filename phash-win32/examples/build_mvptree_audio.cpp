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

  double maxC = 0.0f;
  for (int i=0;i<Nc;i++){
	  if (ptrC[i] > maxC){
          maxC = ptrC[i];
	  }
  }  
  if (ptrC != NULL) 
	  delete [] ptrC;

  return (float)(1000*(1-maxC));

}

int main(int argc, char **argv){
	if (argc < 3){
       printf("not enough input args\n");
       return 1;
	}
 
    const char *dir_name = argv[1];/* name of dir to retrieve image files */
    const char *filename = argv[2];/* name of file to save db */

    float alpha = 2.0f;
    float lvl = 1.0f;

    MVPFile mvpfile;
    ph_mvp_init(&mvpfile);
    mvpfile.leafcapacity = 10;
    mvpfile.pgsize = 1 << 17; /* 131,072 */
    mvpfile.filename = strdup(filename);
    mvpfile.hashdist = distancefunc;
    mvpfile.hash_type = UINT32ARRAY;


    int nbfiles = 0;
    printf("dir name: %s\n", dir_name);
    char **files = ph_readfilenames(dir_name,nbfiles);
    if (!files){
		printf("mem alloc error\n");
	exit(1);
    }
    printf("nbfiles = %d\n", nbfiles);
    DP **hashlist = new DP*[nbfiles];
    if (!hashlist){
	    printf("mem alloc error\n");
	    exit(1);
    }
    
    int count = 0;
    float *sigbuf = new float[1<<20];
    int buflen = (1<<20);
 
    uint32_t *hashspace = new uint32_t[1<<18];
    int hashspacelength = (1<<18);
 
    float *buf;
    uint32_t *hash = hashspace;
    int hashspaceleft = hashspacelength;
    int nbframes;
    for (int i=0;i<nbfiles;i++){
		printf("file[%d]: %s\n", i, files[i]);
        hashlist[count] = ph_malloc_datapoint(mvpfile.hash_type,mvpfile.pathlength);
		if (hashlist[count] == NULL){
			printf("mem alloc error\n");
			continue;
		}
		hashlist[count]->id = files[i];
        buflen = (1<<20);
        buf = ph_readaudio(files[i], 8000, 1, sigbuf, buflen);
		if (buf == NULL){
            printf("unable to get signal\n");
            continue;
		}
        printf("sig buffer length %d\n", buflen);

        hash = ph_audiohash(buf, buflen, hash, hashspaceleft, 8000, nbframes);
		hashlist[count]->hash = hash;
		if (hashlist[count]->hash == NULL){
			printf("unable to get hash\n\n");
			continue;
		}
        printf("hash length %d\n", nbframes);
		hashlist[count]->hash_length = (uint16_t)nbframes;
        
        hash += nbframes;
        hashspaceleft -= nbframes;

        count++;
    }
 
    MVPRetCode errcode = ph_save_mvptree(&mvpfile, hashlist, count);
    if (errcode != PH_SUCCESS){
		printf("unable to save %s, err %d\n", filename, errcode);
		return 1;
    }
    printf("saved files\n");

    
	for (int i=0;i<nbfiles;i++){
        delete [] files[i];
        delete hashlist[i];
	}
    delete [] files;
    delete [] hashlist;
    delete [] sigbuf;
    delete [] hashspace;

    return 0;
}
