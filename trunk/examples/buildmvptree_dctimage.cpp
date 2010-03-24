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
    D Grant Starkweather - dstarkweather@phash.org

*/
#include <stdio.h>
#include <math.h>
#include "pHash.h"


float distancefunc(DP *pa, DP *pb){
    float d = ph_hamming_distance(*((ulong64*)pa->hash), *((ulong64*)pb->hash));
    return d;
}

int main(int argc, char **argv){
    if (argc < 3){
	printf("not enough input args\n");
        printf("usage: %s dirname filename\n", argv[0]);
	return -1;
    }
 
    const char *dir_name = argv[1];/* name of dir to retrieve image files */
    const char *filename = argv[2];/* name of file to save db */

    MVPFile mvpfile;
    mvpfile.branchfactor = 2;
    mvpfile.pathlength = 5;
    mvpfile.leafcapacity = 50;
    mvpfile.pgsize = 8192;
    mvpfile.filename = strdup(filename);
    mvpfile.hashdist = distancefunc;
    mvpfile.hash_type =  UINT64ARRAY;

    int nbfiles = 0;
    printf("dir name: %s\n", dir_name);
    char **files = ph_readfilenames(dir_name,nbfiles);
    if (!files){
	printf("mem alloc error\n");
	return -2;
    }
    printf("nbfiles = %d\n", nbfiles);
    DP **hashlist = (DP**)malloc(nbfiles*sizeof(DP*));
    if (!hashlist){
	printf("mem alloc error\n");
	return -3;
    }

    int count = 0;
    for (int i=0;i<nbfiles;i++){
        hashlist[count] = ph_malloc_datapoint(mvpfile.hash_type);
	if (hashlist[count] == NULL){
	    printf("mem alloc error\n");
	    return -4;
	}
	hashlist[count]->hash = malloc(sizeof(ulong64));
	if (hashlist[count]->hash == NULL){
	    printf("mem alloc error\n");
	    return -5;
	}

        if (ph_dct_imagehash(files[i],*((ulong64*)hashlist[count]->hash)) < 0){
	    printf("unable to get hash\n");
            continue;
	}

        printf("files[%d]: %s hash = %llx\n", i, files[i], *((ulong64*)hashlist[count]->hash));

        hashlist[count]->id = files[i];
	hashlist[count]->hash_length = 1;
        count++;
    }
 
    
    MVPRetCode ret = ph_save_mvptree(&mvpfile, hashlist, count);
    printf("save: ret code %d\n", ret);



    for (int i=0;i<nbfiles;i++){
	free(files[i]);
    }
    free(files);

    for (int i=0;i<nbfiles;i++){
	free(hashlist[i]->hash);
	ph_free_datapoint(hashlist[i]);
    }
    free(hashlist);

    return 0;
}
