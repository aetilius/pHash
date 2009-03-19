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

#include <cstdio>
#include <cstdlib>
#include <dirent.h>
#include <errno.h>
#include <algorithm>
#include <vector>
#include <string.h>
#include "audiophash.h"

#define TRUE 1
#define FALSE 0

using namespace std;

/**
 * struct to hold hash
 **/
struct data_point {
    int nbframes;    //number of frames, i.e. length of hash
    char *name;      //name of file from which hash is computed
    uint32_t *hash;  //pointer to start of hash
};

struct data_point* ph_malloc_dp(){
    
    return (struct data_point*)malloc(sizeof(struct data_point));
    
}
/**
 * purpose: auxiliary function to compare data points by name for sorting
 **/
bool cmp_lt_dp(struct data_point dpa, struct data_point dpb){
    int result = strcmp(dpa.name, dpb.name);
    if (result < 0)
	return TRUE;
    return FALSE;
}
/**
 *  TEST
 * 
 * purpose: test harness for audiohash functions; reads in files from two directories that are
 *          specified on the command line;  same number of files in each directory; each file name
 *          represented in each directory.  The program reads in the files, computes the hash 
 *          and sorts the hash into a vector. Each hash is then compared to the corresponding hash
 *          from the other directory in the sorted order.  These are the "inter" comparisons, where
 *          a confidence vector is computed, the max value being used as a confidence score to
 *          indicate the percentage confidence of a similar file.
 *
 *          The program then uses the first vector to compare hashes from different images.  Each
 *          hash is compared to each successive hash in the list in sorted order.  The max value
 *          of the confidence vector is used as a confidence score to indicate the confidence of 
 *          a match.  These are the "inter" comparisons.
 *
 *          params
 *                  ber threshold  - bit error rate threshold (0.25 to 0.35)
 *                  block_size     - size of block of frames to compare at a time (256)
 *      
 *          
 *         Confidence scores range from 0 to 1.0.  A good threshold for determining similarity
 *         is 0.50.
**/
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
    int N;
    int nbframes;
    float *buf = NULL;
    uint32_t *hash = NULL;
    errno = 0;
    struct dirent *dir_entry;
    struct data_point *dp;
    vector<struct data_point> hashlist;  //vector for first directory hashes
    vector<struct data_point> hashlist2; //vector for second directory hashes
    int nbfiles = 0;
    unsigned int i=0;


    printf("reading files in %s ...\n",dir_name);
    //open first directory
    DIR *dir = opendir(dir_name);
    if (!dir){
	printf("could not open directory\n");
	exit(1);
    }
    char path[100];
    path[0]='\0';
    while ((dir_entry = readdir(dir))!=0) {
	if (strcmp(dir_entry->d_name,".") && strcmp(dir_entry->d_name,"..")){
	    strcat(path,dir_name);
	    strcat(path,"/");
	    strcat(path,dir_entry->d_name);
	    buf = ph_readaudio(path,sr,channels,N);
	    if (!buf)
		continue;
	    hash = ph_audiohash(buf,N,sr,nbframes);
            if (!hash)
                continue;
	    dp = ph_malloc_dp();
	    dp->name = strdup(dir_entry->d_name);
	    dp->hash = hash;
	    dp->nbframes = nbframes;
	    hashlist.push_back(*dp);
	    i++;
	    nbfiles++;
	}
	path[0]='\0';
	errno = 0;
    }

    printf("sorting list ...\n");
    sort(hashlist.begin(),hashlist.end(),cmp_lt_dp);

    printf("hit any key\n");
    getchar();

    printf("reading files in %s ... \n",dir_name2);
    //open first directory
    dir = opendir(dir_name2);
    if (!dir){
	printf("could not open directory\n");
	exit(1);
    }
    i=0;
    path[0]='\0';
    while ((dir_entry = readdir(dir)) != 0) {
	if (strcmp(dir_entry->d_name,".") && strcmp(dir_entry->d_name,"..")){
	    strcat(path,dir_name2);
	    strcat(path,"/");
	    strcat(path,dir_entry->d_name);
	    buf = ph_readaudio(path,sr,channels,N);
	    if (!buf)
		continue;
	    hash = ph_audiohash(buf,N,sr,nbframes);
            if (!hash)
                continue;
	    dp = ph_malloc_dp();
	    dp->name = strdup(dir_entry->d_name);
	    dp->hash = hash;
	    dp->nbframes = nbframes;
	    hashlist2.push_back(*dp);
	    i++;
	    nbfiles++;
	}
	path[0]='\0';
	errno = 0;
    }

    if (errno){
	printf("error reading files\n");
	exit(1);
    }
    
    printf("sorting list ...\n");
    sort(hashlist2.begin(),hashlist2.end(),cmp_lt_dp);

    printf("hit any key\n");
    getchar();

    int Nc;
    double maxC = 0.0;
    double *pC = NULL;

    printf("intra confidence values:\n");
    printf("************************\n");
    for (i=0;i<hashlist.size();i++){
	printf("file1: %s    file2: %s ",hashlist[i].name,hashlist2[i].name);
	pC = ph_audio_distance_ber(hashlist[i].hash, hashlist[i].nbframes, hashlist2[i].hash, hashlist2[i].nbframes, threshold, block_size, Nc);

	maxC = 0.0;
	for (int j=0;j<Nc;j++){
	    if (pC[j] > maxC){
		maxC = pC[j];
	    }
	}
	printf("cs = %3.2f\n",maxC);
    }

    printf("***********************\n");
    printf("hit any key\n");
    getchar();


    printf("inter confidence values\n");
    printf("***********************\n");
    for (i=0;i<hashlist.size();i++){
	for (unsigned int j=i+1;j<hashlist.size();j++){
	    printf("file1:%s file2: %s ",hashlist[i].name,hashlist[j].name);
	    pC = ph_audio_distance_ber(hashlist[i].hash, hashlist[i].nbframes,hashlist[j].hash,hashlist[j].nbframes, threshold, block_size, Nc);
	    maxC = 0.0;
	    for (int k=0;k<Nc;k++){
		if (pC[k]>maxC)
		    maxC = pC[k];
	    }
	    printf("cs=%3.2f\n",maxC);
	}
    }
    printf("done\n");

    return 0;
}
