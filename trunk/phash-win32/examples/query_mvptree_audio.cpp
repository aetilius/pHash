#include <stdio.h>
#include "pHash.h"
#include "audiophash.h"

static int nb_calcs;

float audiohashdistance(DP *dpA, DP *dpB){

    nb_calcs++;

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
    float d = 10*(1-maxC);
    float res = exp(d)-1;
    return res;
}

int main(int argc, char **argv){
	if (argc < 3){
        printf("not enough input args\n");
        exit(1);
	}
   
    const char *dir_name = argv[1];/* name of files in directory of query images */
    const char *filename = argv[2];/* name of file to save db */

    const int sr = 8000;
    const int nbchannels = 1;

    MVPFile mvpfile;
    ph_mvp_init(&mvpfile);
    mvpfile.filename = strdup(filename);
    mvpfile.hashdist = audiohashdistance;
    mvpfile.hash_type = UINT32ARRAY;
   
    int nbfiles = 0;
    printf("using db %s\n", filename);
    printf("using dir %s for query files\n", dir_name);
    char **files = ph_readfilenames(dir_name,nbfiles);
    if (!files){
		printf("mem alloc error\n");
		exit(1);
    }

    printf("nb query files = %d\n", nbfiles);

    DP *query = NULL;
    float radius = 1000.0f;
    if (argc >= 4) radius = atof(argv[3]);
    int knearest = 20;
    if (argc >= 5) knearest = atoi(argv[4]);
    float threshold = 400.0f;
    if (argc >=6) threshold = atof(argv[5]);

    DP **results = (DP**)malloc(knearest * sizeof(DP**));
    if (!results){
		printf("mem alloc error\n");
		exit(1);
    }
    
	int nbfound = 0, count = 0, sum_calcs = 0;
    float *buf;
    int N = 0;
    uint32_t *hash;
    int hashlen = 0;
    printf("******************************\n");
    for (int i=0;i<nbfiles;i++){
		printf("file[%d]: %s\n", i, files[i]);
		buf = ph_readaudio(files[i], sr, nbchannels, NULL, N);
		if (!buf){
			printf("could not read audio\n");
			continue;
		}
		hash = ph_audiohash(buf,N,NULL, 0, sr,hashlen);
		if (!hash){
			printf("could not get hash\n");
			free(buf);
			continue;
		}
        query = ph_malloc_datapoint(mvpfile.hash_type,mvpfile.pathlength);
        if (!query){
			printf("mem alloc error\n");
            free(buf);
			free(hash);
			break;
		}
		query->id = strdup(files[i]);
		query->hash = hash;
		query->hash_length = hashlen;

		count++;
		nb_calcs = 0;
		nbfound = 0;
		MVPRetCode ret = ph_query_mvptree(&mvpfile,query,knearest,radius,threshold, results,&nbfound);
		if (ret != PH_SUCCESS){
			printf("could not complete query\n");
			continue;
		}
		sum_calcs += nb_calcs;

		printf(" %d found, %d distance calcs\n", nbfound,nb_calcs);
		for (int i=0;i<nbfound;i++){
			printf("    %d  %s dist = %f\n", i, results[i]->id, audiohashdistance(results[i], query));
		}
        printf("******************************************\n");
		free(query);
		free(buf);
		free(hash);
    } 
   float ave_calcs = (float)sum_calcs/(float)count;      
   printf("ave calcs/query: %f\n", ave_calcs);
    
    return 0;
}
