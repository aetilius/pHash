#include <stdio.h>
#include <math.h>
#include "pHash.h"
#include "audiophash.h"

static int nb_calcs; /* metric to evaluate efficiency */
                     /* set to 0 before calling build, add, query */  

float distfunc(DP *dpA, DP *dpB){

    nb_calcs++;

    uint32_t *hash1 = (uint32_t*)dpA->hash;
    int N1 = dpA->hash_length;
    uint32_t *hash2 = (uint32_t*)dpB->hash;
    int N2 = dpB->hash_length;

    float threshold = 0.30f;
    int blocksize = 256;
    int Nc=0;
    double *ptrC = ph_audio_distance_ber(hash1, N1, hash2, N2, threshold,blocksize,Nc);

    float maxC = 0;
    for (int i=0;i<Nc;i++){
		if (ptrC[i] > maxC)
			maxC = (float)ptrC[i];
	}
    if (ptrC) free(ptrC);

    double res = 1000*(1 - maxC);
    return res;
}

int main(int argc, char **argv){
	if (argc < 3){
        printf("not enough input args\n");
		printf("usage: progname querydirectory dbname searchradius rnearest threshold\n");
        exit(1);
	}
   
    const char *dir_name = argv[1];/* name of files in directory of query images */
    const char *filename = argv[2];/* name of file to save db */

    const int sr = 8000;        /* sample to convert all files */ 
    const int nbchannels = 1;   /* number channels to convert all files */ 
    const float nbsecs = 45.0f; /* first number seconds to take of query file */ 

    MVPFile mvpfile;            /* mvp configuration */ 
    mvpfile.filename = strdup(filename);
    mvpfile.hashdist = distfunc;
    mvpfile.hash_type = UINT32ARRAY;
   
    /* get list of filenames in directory */ 
    int nbfiles = 0;
    printf("using db %s\n", filename);
    printf("using dir %s for query files\n", dir_name);
    char **files = ph_readfilenames(dir_name,nbfiles);
    if (!files){
		printf("mem alloc error\n");
		exit(1);
    }

    printf("nb query files = %d\n", nbfiles);

    float radius = 700.0f; /* search radius to traverse the tree */ 
    if (argc >= 4) radius = atof(argv[3]);
    int rnearest = 1;      /* max number results to retrieve */ 
    if (argc >= 5) rnearest = atoi(argv[4]);
    float threshold = 300.0; /*distance threshold for inclusion into list */ 
    if (argc >=6) threshold = atof(argv[5]);

    /* get hashes for querys - multithreaded */ 
    DP **querys = ph_audio_hashes(files,nbfiles,sr,nbchannels);
	if (querys == NULL){
        printf("unable to get hashes for querys\n");
        return -1;
	}
    /* check for null entries and move up in list */ 
    int count = 0;   
	for (int i=0;i<nbfiles;i++){
		if (!(querys[i] == NULL  || querys[i]->hash == NULL)){
             querys[count] = querys[i];
             count++;
		}
	}
    DP **results = (DP**)malloc(rnearest*sizeof(DP*));
	if (results == NULL){
        printf("unable to alloc memory\n");
        return -1;
	}

    /* perform querys for each */ 
    int nbfound = 0, sum_calcs = 0;
    for (int i=0;i<count;i++){
		printf("file[%d]: %s hash length %d\n", i, querys[i]->id,querys[i]->hash_length);
		nb_calcs = 0;
		nbfound = 0;
        printf("do query ...\n");
		MVPRetCode ret = ph_query_mvptree(&mvpfile,querys[i],rnearest,radius,threshold, results,&nbfound);
		if (ret != PH_SUCCESS && ret != PH_ERRCAP){
			printf("could not complete query - %d\n",ret);
			continue;
		}
		sum_calcs += nb_calcs;
		printf(" %d found, %d distance calcs\n", nbfound,nb_calcs);
		for (int j=0;j<nbfound;j++){
			printf("    %d  %s dist = %f\n", i, results[j]->id, distfunc(results[j], querys[i]));
		}
        printf("******************************************\n");


		for (int j=0;j<nbfound;j++){
            ph_free_datapoint(results[j]);
		}
    } 
   float ave_calcs = (float)sum_calcs/(float)count;      
   printf("ave calcs/query: %f\n", ave_calcs);
    
   for (int j=0;j<count;j++){
       ph_free_datapoint(querys[j]);
   }
   free(querys);
   free(results);

    return 0;
}
