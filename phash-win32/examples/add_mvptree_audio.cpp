#include <stdio.h>
#include "pHash.h"
#include "audiophash.h"

float audiohashdistance(DP *dpA, DP *dpB){

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
    free(ptrC);
    return (float)maxC;
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
    float *buf = NULL;
    float *sigbuf = (float*)malloc(1<<28);
    int buflen = 1<<26;
    int N;
    uint32_t *hash;
    int nbframes = 0;
    for (int i=0;i<nbfiles;i++){
		printf("file[%d]: %s\n", i, files[i]);
        N = buflen;
		buf = ph_readaudio(files[i], sr, nbchannels, sigbuf, N);
		if (!buf){
			printf("unable to read file\n");
			continue;
		}
		hash = ph_audiohash(buf, N, NULL, 0, sr, nbframes);
		if (!hash){
			printf("unable to get hash\n");
			free(buf);
		}
		printf("nb hashes = %d\n", nbframes);

        hashlist[count] = ph_malloc_datapoint(mvpfile.hash_type,mvpfile.pathlength);
		if (!hashlist[count]){
			printf("mem alloc error\n");
			free(buf);
			free(hash);
            continue;
		}
		hashlist[count]->id = files[i];
		hashlist[count]->hash = (void*)hash;
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
