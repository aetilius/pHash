#include <stdio.h>
#include "dirent.h"
#include "pHash.h"



int main(int argc, char **argv){
 
	if (argc < 3){
         printf("not enough input args\n");
         exit(1);
	}

    const char *dir_name = argv[1];/* name of dir to retrieve image files */
    const char *filename = argv[2];/* name of file to save db */

    MVPFile mvpfile;
    ph_mvp_init(&mvpfile);
    mvpfile.filename = strdup(filename);
    mvpfile.hashdist = hammingdistance;
    mvpfile.hash_type = UINT64ARRAY;


    ulong64 tmphash = 0;
    
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
    for (int i=0;i<nbfiles;i++){
		printf("file[%d]: %s ", i, files[i]);
		if (ph_dct_imagehash(files[i],tmphash) < 0){
			printf("could not get hash for file %s\n",files[i]);
			continue;
		}
		printf("hash = %llx\n", tmphash);

        hashlist[count] = ph_malloc_datapoint(mvpfile.hash_type,mvpfile.pathlength);
        if (!hashlist[count]){
			printf("mem alloc error\n");
			exit(1);
		}
		hashlist[count]->id = files[i];
        void *ptr_hash = malloc(8);
		if (!ptr_hash){
            printf("unable to allocate mem\n");
            exit(1);
        }
        hashlist[count]->hash = ptr_hash;
		ulong64 *ptr = (ulong64*)hashlist[count]->hash;
        *ptr = tmphash;
		hashlist[count]->hash_length = 1;
        count++;
    }

 
    printf("save mvp tree in %s\n", mvpfile.filename);
    if (ph_save_mvptree(&mvpfile, hashlist, count) < 0){
		printf("unable to save %s\n", filename);
		exit(1);
    }

    return 0;
}
