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


#include <stdio.h>
#include <stdlib.h>
#include "pHash.h"

int main(int argc, char **argv){

    if (argc < 3){
	printf("not enough input args\n");
	exit(1);
    }
    const char *file1 = argv[1];
    const char *file2 = argv[2];

    printf("file1: %s\n", file1);
    int nbhashes1 = 0;
    TxtHashPoint *hash1 = ph_texthash(file1,&nbhashes1);
    if (!hash1){
	printf("unable to complete text hash function\n");
	exit(1);
    }
    printf("length %d\n", nbhashes1);

    
    printf("file2: %s\n", file2);
    int nbhashes2 = 0;
    TxtHashPoint *hash2 = ph_texthash(file2,&nbhashes2);
    if (!hash2){
	printf("unable to complete text hash function\n");
	exit(2);
    }

    printf("length %d\n", nbhashes2);
    int count,j;
    TxtMatch *matches = ph_compare_text_hashes(hash1, nbhashes1, hash2, nbhashes2, &count);
    if (!matches){
	printf("unable to complete compare function\n");
	exit(3);
    }
    
    printf(" %d matches\n", count);
    printf(" indxA  indxB  length\n");
    for (j=0;j<count;j++){
	printf(" %jd %jd %d\n", matches[j].first_index, matches[j].second_index,matches[j].length); 
    }

    return 0;
}
