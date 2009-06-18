/*

    pHash, the open source perceptual hash library
    Copyright (C) 2008-2009 Aetilius, Inc.
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

#ifndef _PHASH_H
#define _PHASH_H

#define cimg_debug 0
#define cimg_display 0

#include <limits.h>
#include <math.h>
#include <dirent.h>
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <string.h>
#include <stdint.h>
#include "CImg.h"
#include "config.h"

using namespace cimg_library;
using namespace std;

#define SQRT_TWO 1.4142135623730950488016887242097

#define ROUNDING_FACTOR(x) (((x) >= 0) ? 0.5 : -0.5) 

#if defined( _MSC_VER) || defined(_BORLANDC_)
typedef unsigned _uint64 ulong64;
typedef signed _int64 long64;
#else
typedef unsigned long long ulong64;
typedef signed long long long64;
#endif

const int MaxFileSize = (1<<20); /* 1GB file size limit (for mvp files) */
const off_t HeaderSize = 64;     /* header size for mvp file */

typedef enum ph_mvp_retcode {
    PH_SUCCESS = 0,   /* success */
    PH_ERRPGSIZE,     /* page size error */
    PH_ERRFILE,       /* file operations */
    PH_ERRMAP,        /* mmap'ing error */
    PH_NOSAVEMVP      /* could not save mvp file */
    PH_ERR_ARGLIST,   /* null arg */
    PH_ERR_NODISTFUNC, /* no dist function in mvpfile structure */
    PH_MEMALLOC,       /* mem alloc error - not enough available memory */
    PH_ERR_NTYPE,      /* unrecognized node type */
    PH_RESULTSFULL     /* more results found than can be supported in ret array */
}MVPRetCode;

typedef enum ph_hashtype {
    BYTEARRAY   = 1,          /* refers to bitwidth of the hash value */
    UINT16ARRAY = 2,
    UINT32ARRAY = 4,
    UINT64ARRAY = 8,
}HashType;

typedef struct ph_file_offset {
    off_t offset;
    uint8_t fileno;
} FileIndex;

/* call back function for mvp tree functions - to performa distance calc.'s*/
typedef float (*hash_compareCB)(DP *pointA, DP *pointB);


typedef struct ph_mvp_file {
    const char *filename;   /* name of db to use */
    char *buf;
    off_t file_pos;
    int fd;
    uint8_t nbdbfiles;
    uint8_t branchfactor; /*branch factor of tree, M(=2)*/

    /*length of path to store distances from vantage points in the struct data_point. 
      used when querying or constructing the tree, P(=5) */
    uint8_t pathlength;  

    uint8_t leafcapacity; /*maximum number of data points to a leaf, K(=25) */

    uint8_t isleaf;       /*boolean flag used in query function */

    /*size of page size to store a leaf of the tree structure.
      Must be >= system page size.  Might need to increase above 
      the system page size to fit all the data points in a leaf.
       internal_pgsize is the same thing but for internal nodes. 
       Set to 0 to use the system pg_size */
    off_t leaf_pgsize;
    off_t internal_pgsize;
    HashType hash_type;

    /*callback function to use to calculate the distance between 2 datapoints */
    hash_compareCB hashdist;

} MVPFile ;


/* convenience function to set var's of mvp tree */
void ph_mvp_init(MVPFile *m){
    m->branchfactor = 2;
    m->pathlength = 5;
    m->leafcapacity = 25;
    m->leaf_pgsize = sysconf(_SC_PAGE_SIZE);     /* use host page size */
    m->internal_pgsize = sysconf(_SC_PAGE_SIZE);
    return;
}


/* structure for a single hash */
typedef struct datapoint {
    char *id;
    void *hash;
    float *path;
    uint16_t hash_length;
    uint8_t hash_type;
}DP;

/*! /brief Radon Projection info
 */
typedef struct ph_projections {
    CImg<uint8_t> *R;           //contains projections of image of angled lines through center
    int *nb_pix_perline;        //the head of int array denoting the number of pixels of each line
    int size;                   //the size of nb_pix_perline
}Projections;

/*! /brief feature vector info
 */
typedef struct ph_feature_vector {
    double *features;           //the head of the feature array of double's
    int size;                   //the size of the feature array
}Features;

/*! /brief Digest info
 */
typedef struct ph_digest {
    char *id;                   //hash id
    uint8_t *coeffs;            //the head of the digest integer coefficient array
    int size;                   //the size of the coeff array
} Digest;


/* /brief alloc a single data point
 *  allocates path array, does nto set id or path
 */
DP* ph_malloc_datapoint(int hashtype, int pathlength);

/** /brief free a datapoint and its path
 *
 */
void ph_free_datapoint(DP *dp);

/*! /brief copyright information
 */
const char* ph_about();

/*! /brief radon function
 *  Find radon projections of N lines running through the image center for lines angled 0
 *  to 180 degrees from horizontal.
 *  /param img - CImg src image
 *  /param  N  - int number of angled lines to consider.
 *  /param  projs - (out) Projections struct 
 *  /return int value - less than 0 for error
 */
int ph_radon_projections(const CImg<uint8_t> &img,int N,Projections &projs);

/*! /brief feature vector
 *         compute the feature vector from a radon projection map.
 *  /param  projs - Projections struct
 *  /param  fv    - (out) Features struct
 *  /return int value - less than 0 for error
*/
int ph_feature_vector(const Projections &projs,Features &fv);

/*! /brief dct 
 *  Compute the dct of a given vector
 *  /param R - vector of input series
 *  /param D - (out) the dct of R
 *  /return  int value - less than 0 for error
*/
int ph_dct(const Features &fv, Digest &digest);

/*! /brief cross correlation for 2 series
 *  Compute the cross correlation of two series vectors
 *  /param x - Digest struct
 *  /param y - Digest struct
 *  /param pcc - double value the peak of cross correlation
 *  /param threshold - double value for the threshold value for which 2 images
 *                     are considered the same or different.
 *  /return - int value - 1 (true) for same, 0 (false) for different, < 0 for error
 */

int ph_crosscorr(const Digest &x,const Digest &y,double &pcc, double threshold = 0.90);

/*! /brief image digest
 *  Compute the image digest for an image given the input image
 *  /param img - CImg object representing an input image
 *  /param sigma - double value for the deviation for a gaussian filter function 
 *  /param gamma - double value for gamma correction on the input image
 *  /param digest - (out) Digest struct
 *  /param N      - int value for the number of angles to consider. 
 *  /return       - less than 0 for error
 */
int ph_image_digest(const CImg<uint8_t> &img,double sigma, double gamma,Digest &digest,int N=180);

/*! /brief image digest
 *  Compute the image digest given the file name.
 *  /param file - string value for file name of input image.
 *  /param sigma - double value for the deviation for gaussian filter
 *  /param gamma - double value for gamma correction on the input image.
 *  /param digest - Digest struct
 *  /param N      - int value for number of angles to consider
 */
int ph_image_digest(const char *file, double sigma, double gamma, Digest &digest,int N=180);


/*! /brief compare 2 images
 *  /param imA - CImg object of first image 
 *  /param imB - CImg object of second image
 *  /param pcc   - (out) double value for peak of cross correlation
 *  /param sigma - double value for the deviation of gaussian filter
 *  /param gamma - double value for gamma correction of images
 *  /param N     - int number for the number of angles of radon projections
 *  /param theshold - double value for the threshold
 *  /return int 0 (false) for different images, 1 (true) for same image, less than 0 for error
 */
int ph_compare_images(const CImg<uint8_t> &imA,const CImg<uint8_t> &imB,double &pcc, double sigma = 3.5, double gamma = 1.0,int N=180,double threshold=0.90);

/*! /brief compare 2 images
 *  Compare 2 images given the file names
 *  /param file1 - char string of first image file
 *  /param file2 - char string of second image file
 *  /param pcc   - (out) double value for peak of cross correlation
 *  /param sigma - double value for deviation of gaussian filter
 *  /param gamma - double value for gamma correction of images
 *  /param N     - int number for number of angles
 *  /return int 0 (false) for different image, 1 (true) for same images, less than 0 for error
 */
int ph_compare_images(const char *file1, const char *file2,double &pcc, double sigma = 3.5, double gamma=1.0, int N=180,double threshold=0.90);

/*! /brief return dct matrix, C
 *  Return DCT matrix of sqare size, N
 *  /param N - int denoting the size of the square matrix to create.
 *  /return CImg<double> size NxN containing the dct matrix
 */
CImg<float>* ph_dct_matrix(const int N);

/*! /brief compute dct robust image hash
 *  /param file string variable for name of file
 *  /param hash of type ulong64 (must be 64-bit variable)
 *  /return int value - -1 for failure, 1 for success
 */
int ph_dct_imagehash(const char* file,ulong64 &hash);

/* ! /brief hamming distance
 *   Compute the hamming distance between two 64-bit data types
 *   /param hash1 of type ulong64 denoting an image hash
 *   /param hash2 of type ulong64 denoting an image hash
 *   /return int value for the hamming distance, -1 for error
 */
/* ! /brief dct video robust hash
 *   Compute video hash based on the dct of normalized video 32x32x64 cube
 *   /param file name of file
 *   /param hash ulong64 value for hash value
 *   /return int value - less than 0 for error
 */
int ph_dct_videohash(const char* file,ulong64 &hash);

int ph_hamming_distance(const ulong64 hash1,const ulong64 hash2);

int ph_rash_videodigest(const char* file,CImg<uint8_t> *p_videodigest);


/** /brief create a list of datapoint's directly from a directory of image files
 *  /param dirname - path and name of directory containg all image file names
 *  /param capacity - int value for upper limit on number of hashes
 *  /param count - number of hashes created (out param)
 *  /return pointer to a list of DP pointers (NULL for error)
 */

DP** ph_read_imagehashes(const char *dirname,int capacity, int &count);

/** /brief get all the filenames in specified directory
 *  /param dirname - string value for path and filename
 *  /param cap - int value for upper limit to number of files
 *  /param count - int value for number of file names returned
 *  /return array of pointers to string file names (NULL for error)
 **/
char** ph_readfilenames(const char *dirname,int cap,int &count);


DP* ph_read_datapoint(MVPFile *m);

off_t ph_save_datapoint(DP *new_dp, MVPFile *m);

MVPFile* _ph_map_mvpfile(uint8_t filenumber, off_t offset, MVPFile *m);

void _ph_unmap_mvpfile(uint8_t filenumber, off_t orig_pos, MVPFile *m, MVPFile *m2);

float hammingdistance(DP *pntA, DP *pntB);

MVPRetCode ph_query_mvptree(MVPFile *m, DP *query, int knearest, float radius,
			    DP **results, int *count, int level);

MVPRetCode ph_query_mvptree(MVPFile *m, DP *query, int knearest, float radius,
			    DP **results, int *count);

FileIndex* ph_save_mvptree(MVPFile *m, DP **points, int nbpoints, int saveall_flag, int level);

MVPRetCode ph_save_mvptree(MVPFile *m, DP **points, int nbpoints);

MVPRetCode ph_add_mvptree(MVPFile *m, DP *new_dp, int level);

int ph_add_mvptree(MVPFile *m, DP **points, int nbpoints);

#endif
