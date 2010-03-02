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

  double res = 1000*(1.0 - maxC);
  return (float)res;

}

int main(int argc, char **argv){
	if (argc < 3){
       printf("not enough input args\n");
	   printf("usage: progname directory dbname\n");
       return 1;
	}
    const char *dir_name = argv[1];/* name of dir to retrieve image files */
    const char *filename = argv[2];/* name of indexing db, e.g. 'audiodb'  */

    const int sr = 8000;            /* convert to sr */
    const int nbchannels = 1;       /* convert to number of channels */
    const float nbsecs = 45.0f;     /* number secs of audio to read from each file */

    MVPFile mvpfile;                /* mvp tree indexing configuration */ 
    mvpfile.branchfactor = 2;       /* number of branches for each node */ 
    mvpfile.pathlength = 5;         /* length of path to remember distances from each level vantage points */
    mvpfile.leafcapacity = 44;      /* number of datapoints for each leaf */ 
	mvpfile.pgsize = 1 << 19;       /* 2^29=524,288 Must be able to fit leafcapacity*sizeof(datapoint) */
    mvpfile.filename = strdup(filename);
    mvpfile.hashdist = distancefunc;    /*distance function to use in the tree */ 
    mvpfile.hash_type = UINT32ARRAY;    /* bit width of each haash element */ 

    /* read in file names from directory */ 
    int nbfiles = 0;
    printf("dir name: %s\n", dir_name);
    char **files = ph_readfilenames(dir_name,nbfiles);
    if (!files){
		printf("mem alloc error\n");
		return 1;
    }
   
   /* retrieve audio hashes - multithreaded */ 
   int count = 0;
   DP** hashlist = ph_audio_hashes(files, nbfiles, sr, nbchannels);  
   if (hashlist == NULL){
       printf("unable to get audio hashes\n");
       return 1;
   }
   /* if any datapoints are null, move the list up */
   for (int i=0;i<nbfiles;i++){
	   if ( !((hashlist[i] == NULL) || (hashlist[i]->hash == NULL))){
          hashlist[count] = hashlist[i];
          count++;
	   }
   }

    /* save tree */ 
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
    return 0;
}
