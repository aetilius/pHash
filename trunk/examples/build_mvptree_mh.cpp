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
    if ((!pa)||(!pb)){
	printf("dp is null\n");
	return -1.0;
    }
    uint8_t *hashA = (uint8_t*)pa->hash;
    uint8_t *hashB = (uint8_t*)pb->hash;
    float d = 10*ph_hammingdistance2(hashA,pa->hash_length,hashB,pb->hash_length);
    float res = exp(d)-1;
    return res;
}

int main(int argc, char **argv){
    if (argc < 3){
	printf("not enough input args\n");
	return 1;
    }
 
    const char *dir_name = argv[1];/* name of dir to retrieve image files */
    const char *filename = argv[2];/* name of file to save db */

    float alpha = 2.0f;
    float lvl = 1.0f;

    MVPFile mvpfile;
    mvpfile.branchfactor = 2;
    mvpfile.pathlength = 5;
    mvpfile.filename = strdup(filename);
    mvpfile.hashdist = distancefunc;
    mvpfile.hash_type = BYTEARRAY;
    mvpfile.pgsize = 8192;
    mvpfile.leafcapacity = 45;    

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
    int hashlength;
    int count = 0;
    for (int i=0;i<nbfiles;i++){
	printf("file[%d]: %s\n", i, files[i]);
        hashlist[count] = ph_malloc_datapoint(mvpfile.hash_type);
	if (hashlist[count] == NULL){
	    printf("mem alloc error\n");
	    exit(1);
	}
	hashlist[count]->id = files[i];
	hashlist[count]->hash = ph_mh_imagehash(files[i],hashlength, alpha, lvl);
	if (hashlist[count]->hash == NULL){
	    printf("unable to get hash\n");
	    exit(1);
	}
	printf("len %d\n", hashlength);
	hashlist[count]->hash_length = hashlength;
        count++;
    }
 
    
    MVPRetCode ret = ph_save_mvptree(&mvpfile, hashlist, count);
    printf("save: ret code %d", ret);

    return 0;
}
