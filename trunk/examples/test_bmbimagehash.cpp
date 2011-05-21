#include <stdio.h>
#include "pHash.h"

using namespace std;

int main(int argc, char **argv)
{
	if(argc != 2) 
	{
		printf("Error: Invalid command line arguments\n");
		exit(1);
	}

	
	char *fileA = argv[1];
	char *fileB = argv[2];

	printf("--- file 1 --- \n");
	
	BinHash *hash_a;
	bmb_imagehash(fileA,1,&hash_a);

	printf("--- file 2 --- \n");
	BinHash *hash_b;
	bmb_imagehash(fileB,1,&hash_b);

	double dist = ph_hammingdistance2(hash_a->hash, hash_a->bytelength, hash_b->hash, hash_b->bytelength);
	printf("distance = %g\n", dist);

    ph_bmb_free(hash_a);
    ph_bmb_free(hash_b);

    return 0;
}
