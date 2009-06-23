#include <stdio.h>
#include <dirent.h>
#include "pHash.h"


int main(int argc, char **argv){
 
    const char *filename = argv[1];/* name of file to save db */

    MVPFile mvpfile;
    ph_mvp_init(&mvpfile);
    mvpfile.filename = strdup(filename);
    mvpfile.hashdist = hammingdistance;
    mvpfile.hash_type = UINT64ARRAY;


    ulong64 tmphash = 0;

    char c = 'y';
    char queryfile[256];
    float radius = 21.0;
    const int knearest = 20;
    DP *results[knearest];
    int nbfound;

    DP *query = ph_malloc_datapoint(mvpfile.hash_type, mvpfile.pathlength);
    do {
	nbfound = 0;
	printf("enter file:");
	fgets(queryfile, 100, stdin);
	queryfile[strlen(queryfile)-1] = '\0';
	if (ph_dct_imagehash(queryfile, tmphash) < 0)
	    continue;
	query->id = strdup(queryfile);
	query->hash = malloc(UINT64ARRAY);
	ulong64 *ptrHash = (ulong64*)query->hash;
	ptrHash[0] = tmphash; 

	printf("search for: %s\n", queryfile);
	int res;
	if ((res = ph_query_mvptree(&mvpfile, query, knearest, radius, results, &nbfound)) != 0){
	    printf("could not complete query, ret result %d\n",res);
	    continue;
	}
	printf(" %d files found\n", nbfound);
	for (int i=0;i<nbfound;i++){
	    ulong64 *hash = (ulong64*)results[i]->hash;
	    printf(" %d %s %llx\n", i, results[i]->id, *hash);
	}
	printf("again? y/n:");
	c = fgetc(stdin);
	getchar();
    } while (c != 'n');


    return 0;
}
