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
#include <errno.h>
#include "pHash.h"
#include "audiophash.h"

#define TRUE 1
#define FALSE 0



int sort_names(char **names, int L1){

    for (int i=0;i<L1;i++){
	int min = i;
	for (int j=i+1;j<L1;j++){
	    if (strcmp(names[j], names[min]) <= 0)
		min = j;
	}
	if (i != min){
	    char *swap = names[i];
	    names[i] = names[min];
	    names[min] = swap;
	}

    }
    return 0;
}

int main(int argc, char **argv){

    if (argc < 3){
	printf("no args");
	exit(1);
    }
    const char *dir_name  = argv[1];     //first directory
    const char *dir_name2 = argv[2];     //second directory
    const float threshold = 0.30;        //ber threshold (0.25-0.35)
    const int block_size = 256;          //number of frames to compare at a time
    const int sr = 8000;                 //sample rate to convert the stream
    const int channels = 1;              //number of channels to convert stream

    int nbfiles1 = 0;
    char **files1 = ph_readfilenames(dir_name, nbfiles1);
    sort_names(files1, nbfiles1);
    int nbfiles2 = 0;
    char **files2 = ph_readfilenames(dir_name2, nbfiles2);
    sort_names(files2, nbfiles2);

    printf("dir: %s %d\n", dir_name, nbfiles1);
    printf("dir: %s %d\n", dir_name2, nbfiles2);
    if (nbfiles1 != nbfiles2){
	printf("directories do not have same number of files\n");
	exit(1);
    }

    uint32_t **hashes = (uint32_t**)malloc(nbfiles1*sizeof(uint32_t*));
    int *lens = (int*)malloc(nbfiles1*sizeof(int));
    float *buf;
    int buflen;
    uint32_t *hash1, *hash2;
    int hashlen1, hashlen2;
    double *cs;
    int Nc;
    int index, i, j;
    
    printf("intra distances\n");
    printf("***************\n");
    for (index=0;index<nbfiles1;index++){
	printf("file1: %s\n", files1[index]);
	buf = ph_readaudio(files1[index], sr, channels, NULL, buflen);
	if (!buf){
	    printf("unable to read audio\n");
	    continue;
	}
	hash1 = ph_audiohash(buf, buflen, sr, hashlen1);
	if (!hash1){
	    printf("unable to get hash\n");
	    continue;
	}
	hashes[index] = hash1;
	lens[index] = hashlen1;
	free(buf);

	printf("file2: %s\n", files2[index]);
	buf = ph_readaudio(files2[index], sr, channels, NULL, buflen);
	if (!buf) {
	    printf("unable to get audio\n");
	    continue;
	}
        hash2 = ph_audiohash(buf, buflen, sr, hashlen2);
	if (!hash2) {
	    printf("unable to get hash\n");
	    continue;
	}

	cs = ph_audio_distance_ber(hash1, hashlen1, hash2, hashlen2, threshold, block_size, Nc);
	if (!cs){
	    printf("unable to calc distance\n");
	    continue;
	}

	double max_cs = 0.0;
	for (i=0;i<Nc;i++){
	    if (cs[i] > max_cs){
		max_cs = cs[i];
	    }
	}
	printf("max cs %f\n\n", max_cs);

	free(hash2);
	free(buf);
	free(cs);
    }
    
    printf("pause - hit any key\n\n");
    getchar();
    printf("inter dists\n");
    printf("***********\n");
    for (i=0;i<nbfiles1;i++){
	printf("file1: %s\n", files1[i]);
	for (j=i+1;j<nbfiles1;j++){
	    printf("    file2: %s\n", files1[j]);
	    cs=ph_audio_distance_ber(hashes[i],lens[i],hashes[j],lens[j],threshold,block_size,Nc);
	    double max_f = 0.0;
	    for (index=0;index<Nc;index++){
		if (cs[index] > max_f)
		    max_f = cs[index];
	    }
	    printf("    cs = %f\n", max_f);
	}
    }



    return 0;
}
