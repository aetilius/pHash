#include <stdio.h>
#include <math.h>
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

    float maxC = 0;
    for (int i=0;i<Nc;i++){
	if (ptrC[i] > maxC)
	    maxC = (float)ptrC[i];
    }
    float d = 100*(1-maxC);
    float res = exp(d)-1;
    printf("dist %f\n", res);
    return res;
}

int main(int argc, char **argv){
	if (argc < 3){
        printf("not enough input args\n");
        exit(1);
	}
    const char *dir_name = argv[1];/* name of dir to retrieve image files */
    const char *filename = argv[2];/* name of file to save db */

    MVPFile mvpfile;
    mvpfile.branchfactor = 2;
    mvpfile.leafcapacity = 10;
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
    DP **hashlist = (DP**)malloc(nbfiles*sizeof(DP*));
    if (!hashlist){
		printf("mem alloc error\n");
		exit(1);
    }

    int count = 0;
    float *buf = NULL;
    float *sigbuf = (float*)malloc((1<<26)*sizeof(float));
    int buflen = (1 << 27)/sizeof(float);;
	if (sigbuf == NULL) { 
		printf("mem alloc err\n"); exit(1);
	}
    printf("signal buffer %p for %d samples\n", sigbuf, buflen);
    int nbsamples;
    uint32_t *hash;
    int hashlen;
    for (int i=0;i<nbfiles;i++){
		printf("file[%d]: %s\n", i, files[i]);
        nbsamples = buflen;
		buf = ph_readaudio(files[i], sr, nbchannels, sigbuf, nbsamples);
		printf("audio in buf: %p for %d samples\n", buf,nbsamples);
		if (!buf){
			printf("unable to read file: %s\n", files[i]);
            continue;
		}
        
        hash = ph_audiohash(buf,nbsamples,NULL,0,sr,hashlen);
		if (!hash){
            printf("unable to get hash\n");
            continue;
		} 

        hashlist[i] = ph_malloc_datapoint(mvpfile.hash_type,mvpfile.pathlength);
        hashlist[i]->id = files[i];
        hashlist[i]->hash = hash;
        hashlist[i]->hash_length = hashlen;

	}
    
   MVPRetCode retcode = ph_save_mvptree(&mvpfile, hashlist, nbfiles);
   if (retcode != PH_SUCCESS){
       printf("unable to save properly, error code %d\n", retcode);
       return 0;
   } 
   printf("saved: %d\n",retcode);
   return 0;
}
