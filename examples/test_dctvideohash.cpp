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
    David Starkweather - dstarkweather@phash.org

*/


#include "stdio.h"
#include "pHash.h"
#include "CImg.h"
#include "cimgffmpeg.h"

int main(int argc, char **argv){
    if (argc < 2){
	printf("not enough input arguments\n");
	return -1;
    }
    const char *file1 = argv[1];
    const char *file2t = argv[2];

    char *file2 = strdup(file2t);

    if (!file1 || !file2){
	printf("not enough arguments\n");
	exit(1);
    }

    ulong64 *hash1, *hash2;
    int L1, L2;
    double sim;

    printf("file1=%s\n", file1);
    hash1 = ph_dct_videohash(file1, L1);
    if (!hash1){
	exit(2);
    }
    printf("length %d\n", L1);
   for (int i=0;i<L1;i++)
       printf("hash1[%d]=%llx\n", i, hash1[i]);

    do {

	printf("file=%s\n", file2);
	hash2 = ph_dct_videohash(file2, L2);
	if (!hash2){
	    printf("hash 2 is null\n");
	    free(hash1);
	    exit(3);
	}

	printf("length %d\n", L2);

	sim = ph_dct_videohash_dist(hash1, L1, hash2, L2, 21);
	printf("similarity %f\n", sim);

	free(hash2);
	hash2 = NULL;

	file2 = fgets(file2, 80, stdin);
	file2[strlen(file2)-1] = '\0';
    } while (strcmp(file2,"exit"));

    free(hash1);
    free(hash2);
    hash1 = NULL;
    hash2 = NULL;
	free(file2);       
    printf("done\n");
    return 0;
}
 
