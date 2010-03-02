#include <stdio.h>
#include <math.h>
#include "pHash.h"
#include "audiophash.h"

float distfunc(DP *dpA, DP *dpB){

    uint32_t *hash1 = (uint32_t*)dpA->hash;
    int N1 = dpA->hash_length;
    uint32_t *hash2 = (uint32_t*)dpB->hash;
    int N2 = dpB->hash_length;

    float threshold = 0.30f;
    int blocksize = 256;
    int Nc=0;
    double *ptrC = ph_audio_distance_ber(hash1, N1, hash2, N2, threshold,blocksize,Nc);

    double maxC = 0.0f;
    for (int i=0;i<Nc;i++){
	if (ptrC[i] > maxC)
	    maxC = ptrC[i];
    }
    if (ptrC) free(ptrC);
    double res = 1000*(1-maxC);
    return (float)res;
}
 
int main(int argc, char **argv){
	if (argc < 3){
        printf("not enough input args\n");
        return -1;
	}
    const char *dir_name = argv[1];/* name of dir to retrieve image files */
    const char *filename = argv[2];/* name of file to save db */

    const int sr = 8000;          /* sample rate to convert all files */
    const int nbchannels = 1;     /* number of channels to convert all files*/ 
    const float nbsecs = 45.0f;   /* first number seconds of each file to read */ 

    MVPFile mvpfile;              /* mvp tree configuration settings */ 
    mvpfile.filename = strdup(filename);
    mvpfile.hashdist = distfunc;
    mvpfile.hash_type = UINT32ARRAY;

    int nbfiles = 0;
    printf("dir name: %s\n", dir_name);
    char **files = ph_readfilenames(dir_name,nbfiles);
    if (!files){
		printf("mem alloc error\n");
		return -1;
    }
    printf("number of files %d\n", nbfiles);

   /* get audio hashes - multithreaded */ 
   DP **hashlist = ph_audio_hashes(files,nbfiles,sr,nbchannels, 0, nbsecs);
   if (hashlist == NULL){
      printf("could not get audiohashes\n");
      return -1;
   }

   /* if any hashes are null, move list up */ 
   int count = 0;
   for (int i=0;i<nbfiles;i++){
	   if (!(hashlist[i]==NULL || hashlist[i]->hash == NULL)){
          hashlist[count] = hashlist[i];
          count++;
	   }
   }

    /* add each hash to list */ 
    printf("add %d files to file %s\n", count, filename);
    int nbsaved;
    MVPRetCode retcode = ph_add_mvptree(&mvpfile, hashlist, count,nbsaved);
    if (retcode != PH_SUCCESS){
		printf("unable to add points to %s - %d\n", filename,retcode);
    }
    printf("save %d files out of %d\n", nbsaved, count);

    return 0;
}
