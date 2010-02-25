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
    double d = 10*(1-maxC);
    float res = (float)(exp(d)-1.0);

    return res;
}
 
int main(int argc, char **argv){
	if (argc < 3){
        printf("not enough input args\n");
        exit(1);
	}
    const char *dir_name = argv[1];/* name of dir to retrieve image files */
    const char *filename = argv[2];/* name of file to save db */

    const int sr = 8000;
    const int nbchannels = 1;
    const float nbsecs = 45.0f;

    MVPFile mvpfile;
    mvpfile.filename = strdup(filename);
    mvpfile.hashdist = distfunc;
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
    float *buf = NULL;
    float *sigbuf = (float*)malloc(1<<21);
    int buflen = (1<<21)/sizeof(float);
    int N;

    uint32_t *hashes = (uint32_t*)malloc(1<<21);
    int hashspacelength = (1<<21)/sizeof(uint32_t);
    int spaceleft = hashspacelength;
    uint32_t *hash = hashes;
    uint32_t *hash1 = NULL;
    int nbframes = 0;

    for (int i=0;i<nbfiles;i++){
		printf("file[%d]: %s\n", i, files[i]);
        N = buflen;
		buf = ph_readaudio(files[i], sr, nbchannels, sigbuf, N, nbsecs);
		if (!buf){
			printf("unable to read file\n");
			continue;
		}
		hash1 = ph_audiohash(buf, N, hash, spaceleft, sr, nbframes);
		if (!hash1){
			printf("unable to get hash\n");
			free(buf);
		}
        hash += nbframes;
        spaceleft -= nbframes;

        hashlist[count] = ph_malloc_datapoint(mvpfile.hash_type);
		if (!hashlist[count]){
			printf("mem alloc error\n");
            break;
		}
		hashlist[count]->id = files[i];
		hashlist[count]->hash = hash1;
		hashlist[count]->hash_length = nbframes;
        count++;
    }

    printf("add %d files to file %s\n", count, filename);
    int nbsaved;
    MVPRetCode retcode = ph_add_mvptree(&mvpfile, hashlist, count,nbsaved);
    if (retcode != PH_SUCCESS){
		printf("unable to add points to %s - %d\n", filename,retcode);
    }
    printf("save %d files out of %d\n", nbsaved, count);

    return 0;
}
