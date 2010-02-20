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
#include <string.h>
#include <vector>
#include "phash.h"
#include "audiophash.h"

#define TRUE 1
#define FALSE 0

using namespace std;

typedef struct hashpoint {
     char *filename;
     uint32_t *hash;
     int length;
} HP ;

int main(int argc, char **argv){
    const char *msg = ph_about();
    printf(msg);
    printf("\n\n");

    if (argc < 3){
		printf("no args");
	    exit(1);
    }
    const char *dir_name  = argv[1];     //first directory
    const char *dir_name2 = argv[2];     //second directory
    const float threshold = 0.30;        //ber threshold (0.25-0.35)
    const int bs = 256;                 //number of frames to compare at a time
    const int sr = 8000;                 //sample rate to convert the stream
    const int channels = 1;              //number of channels to convert stream

    int N1;
    char **files1 = ph_readfilenames(dir_name, N1);
    int N2;
    char **files2 = ph_readfilenames(dir_name2, N2);
	printf("dir1: %s %d\n", dir_name, N1);
	printf("dir2: %s %d\n", dir_name2, N2);

	if (N1 != N2){
        printf("unequal number files in both directories\n");
        return -1;
	}
    vector<HP> HashList1(N1);

    float *tmpbuf;
    float *buf = (float*)malloc(1<<27);
    int buflen = (1<<27)/sizeof(float);
	if (!buf){
         printf("unable to alloc buffer storage\n");
         return -1;
	}
    

    int hashN = (1 << 26)/sizeof(uint32_t);
    uint32_t *hashes = (uint32_t*)malloc(1 << 26);
	if (!hashes){
        printf("unable to alloc storage for hashes\n");
        return -1;
	}
    int hash_index = 0;
    uint32_t *hash1, *hash2;
    int buflen1, buflen2;
    int nbframes1, nbframes2;
    double *dist;
    int Nc;
    printf("intra distances\n");
	for (int i=0;i<N1;i++){
         printf("  files1[%d] = %s\n", i, files1[i]);
         buflen1 = buflen;
         tmpbuf = ph_readaudio(files1[i], sr, channels, buf, buflen1);
         if (!tmpbuf) continue;
         buf  = tmpbuf;
		 for (int j=0;j<buflen1;j++){
             printf("buf[%d]=%f\n", j, buf[j]);
		 }
        printf("hit any key\n");
        getchar();  
	}
    printf("hit any key to continue\n");
    getchar();
    printf("inter distances\n");
	for (int i=0;i<(int)HashList1.size();i++){
		for (int j=i+1;j<(int)HashList1.size();j++){
			printf("file: %s\n", HashList1[i].filename);
			printf("file: %s\n", HashList1[j].filename);
			dist = ph_audio_distance_ber(HashList1[i].hash, HashList1[i].length, HashList1[j].hash, HashList1[j].length, threshold, bs, Nc);               
            double maxC = 0.0;
			for (int k=0;k<Nc;k++){
                if (dist[k] > maxC)
                    maxC = dist[k];
			}
            printf("cs = %f\n", maxC);
          
            free(dist);

		}
	}



    free(files1);
    free(files2);
    free(hashes);
    free(buf);

    return 0;
}
