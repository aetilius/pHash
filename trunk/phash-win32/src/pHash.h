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

#ifdef WIN32
#define snprintf sprintf_s
#include <share.h>
#endif

#include "CImg.h"
#include "stdint.h"

using namespace cimg_library;
using namespace std;

#define SQRT_TWO 1.4142135623730950488016887242097

#ifndef ULLONG_MAX
#define ULLONG_MAX 18446744073709551615ULL
#endif

#if defined( _MSC_VER) || defined(_BORLANDC_)
#define PACKAGE_STRING "pHash 0.50"
typedef unsigned __int64 ulong64;
typedef signed __int64 long64;
#else
typedef unsigned long long ulong64;
typedef signed long long long64;
#endif

const int MaxFileSize = (1<<30); /* 1GB file size limit (for mvp files) */
const off_t HeaderSize = 64;     /* header size for mvp file */


__declspec(dllexport) const char *mvptag = "pHashMVPfile2009";

typedef enum ph_mvp_retcode {
    PH_SUCCESS = 0,   /* success */
    PH_ERRPGSIZE,     /* page size error */
    PH_ERRFILE,       /* file operations */
    PH_ERRMMAP,        /* mmap'ing error */
    PH_ERRMSYNC,       /* msync error */
    PH_ERRTRUNC,       /* error truncating file */
    PH_ERRSAVEMVP,      /* could not save mvp file */
    PH_ERRARG,   /* null arg */
    PH_ERRMEM,       /* mem alloc error - not enough available memory */
    PH_ERRNTYPE,      /* unrecognized node type */
    PH_ERRCAP,     /* more results found than can be supported in ret array */
    PH_ERRFILETYPE,  /*unrecognized file type  */
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


/* structure for a single hash */
typedef struct ph_datapoint {
    char *id;
    void *hash;
    float *path;
    uint16_t hash_length;
    uint8_t hash_type;
}DP;

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

    /*size of page size to store a leaf of the tree structure.
      Must be >= system page size.  Might need to increase above 
      the system page size to fit all the data points in a leaf.
    */
    off_t pgsize;
    HashType hash_type;

    /*callback function to use to calculate the distance between 2 datapoints */
    hash_compareCB hashdist;

} MVPFile ;


/* convenience function to set var's of mvp tree */
__declspec(dllexport)
void ph_mvp_init(MVPFile *m);

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


/* variables for textual hash */
const int KgramLength = 50;
const int WindowLength = 100;
const int delta = 1;

typedef struct ph_hash_point {
    ulong64 hash;
    off_t index; /*pos of hash in orig file */
} TxtHashPoint;

typedef struct ph_match{
    off_t first_index; /* offset into first file */
    off_t second_index; /* offset into second file */
    uint32_t length;    /*length of match between 2 files */
} TxtMatch;

/* /brief alloc a single data point
 *  allocates path array, does nto set id or path
 */
__declspec(dllexport)
DP* ph_malloc_datapoint(int hashtype, int pathlength);

/** /brief free a datapoint and its path
 *
 */
__declspec(dllexport)
void ph_free_datapoint(DP *dp);

/*! /brief copyright information
 */
__declspec(dllexport) const char* ph_about();

/*! /brief radon function
 *  Find radon projections of N lines running through the image center for lines angled 0
 *  to 180 degrees from horizontal.
 *  /param img - CImg src image
 *  /param  N  - int number of angled lines to consider.
 *  /param  projs - (out) Projections struct 
 *  /return int value - less than 0 for error
 */
__declspec(dllexport)
int ph_radon_projections(const CImg<uint8_t> &img,int N,Projections &projs);

/*! /brief feature vector
 *         compute the feature vector from a radon projection map.
 *  /param  projs - Projections struct
 *  /param  fv    - (out) Features struct
 *  /return int value - less than 0 for error
*/
__declspec(dllexport)
int ph_feature_vector(const Projections &projs,Features &fv);

/*! /brief dct 
 *  Compute the dct of a given vector
 *  /param R - vector of input series
 *  /param D - (out) the dct of R
 *  /return  int value - less than 0 for error
*/
__declspec(dllexport)
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
__declspec(dllexport)
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
__declspec(dllexport)
int ph_image_digest(const CImg<uint8_t> &img,double sigma, double gamma,Digest &digest,int N=180);

/*! /brief image digest
 *  Compute the image digest given the file name.
 *  /param file - string value for file name of input image.
 *  /param sigma - double value for the deviation for gaussian filter
 *  /param gamma - double value for gamma correction on the input image.
 *  /param digest - Digest struct
 *  /param N      - int value for number of angles to consider
 */
__declspec(dllexport)
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
__declspec(dllexport)
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
__declspec(dllexport)
int ph_compare_images(const char *file1, const char *file2,double &pcc, double sigma = 3.5, double gamma=1.0, int N=180,double threshold=0.90);

/*! /brief return dct matrix, C
 *  Return DCT matrix of sqare size, N
 *  /param N - int denoting the size of the square matrix to create.
 *  /return CImg<double> size NxN containing the dct matrix
 */
__declspec(dllexport)
CImg<float>* ph_dct_matrix(const int N);

/*! /brief compute dct robust image hash
 *  /param file string variable for name of file
 *  /param hash of type ulong64 (must be 64-bit variable)
 *  /return int value - -1 for failure, 1 for success
 */
__declspec(dllexport)
int ph_dct_imagehash(const char* file,ulong64 &hash);

/* !calculate a hash for video file
*  /param string filename 
*  /param int Length (out) of hash returned
*  /return ulong64* array of hash values 
*/
__declspec(dllexport)
ulong64* ph_dct_videohash(const char *filename, int &Lenght);

/* !distance function for distance between two video hashes
*  /param hashA - ptr to hashA
*  /param N1 length of hashA
*  /param hashB - ptr to hashB
*  /param N2 length of hashB
*  /param threshold hamming distance threshold to use between each hash value
*  /param double value representing distance, -1 for error
*/
__declspec(dllexport)
double ph_dct_videohash_dist(ulong64 *hashA, int N1, ulong64 *hashB, int N2, int threshold=21);

/* ! /brief dct video robust hash
 *   Compute video hash based on the dct of normalized video 32x32x64 cube
 *   /param file name of file
 *   /param hash ulong64 value for hash value
 *   /return int value - less than 0 for error
 */
__declspec(dllexport)
int ph_hamming_distance(const ulong64 hash1,const ulong64 hash2);

/** /brief create a list of datapoint's directly from a directory of image files
 *  /param dirname - path and name of directory containg all image file names
 *  /param capacity - int value for upper limit on number of hashes
 *  /param count - number of hashes created (out param)
 *  /return pointer to a list of DP pointers (NULL for error)
 */
__declspec(dllexport)
DP** ph_read_imagehashes(const char *dirname,int capacity, int &count);

/** /brief get all the filenames in specified directory
 *  /param dirname - string value for path and filename
 *  /param cap - int value for upper limit to number of files
 *  /param count - int value for number of file names returned
 *  /return array of pointers to string file names (NULL for error)
 **/
__declspec(dllexport)
char** ph_readfilenames(const char *dirname,int &count);

/** /brief read a datapoint (aux function)
 *   /param m current MVPFile struct containing state information
 *   /return DP* the read datapoint struct
 **/
__declspec(dllexport)
DP* ph_read_datapoint(MVPFile *m);

/** /brief get size of a datapoint in bytes (aux. function)
 *  /param m MVPFile struct 
 *  /param dp DP struct
 *  /return int number of bytes the datapoint will use in the file.
 **/
__declspec(dllexport)
int ph_sizeof_dp(DP *dp,MVPFile *m);

/**  /brief save datapoint to file (aux. function)
 *   /param new_dp - DP struct of dp to be saved.
 *   /param m - MVPFil
 *   /return off_t file offset of newly written dp.
 **/
__declspec(dllexport)
off_t ph_save_datapoint(DP *new_dp, MVPFile *m);

/** /brief mmap memory to filenumber/offset
 *  /param filenumber - uint8_t number of file to map
 *  /param offset - off_t offset into new file
 *  /param m - MVPFile
 *  /return MVPFile - ptr to new struct containing the mmap info
 **/
__declspec(dllexport)
MVPFile* _ph_map_mvpfile(uint8_t filenumber, off_t offset, MVPFile *m);

/** /brief unmap/map from m2 to m
 *  /param filenumber - uint8_t filenumber of m2
 *  /param orig_pos   = off_t offset into original file in m.
 *  /return void
 **/
__declspec(dllexport)
void _ph_unmap_mvpfile(uint8_t filenumber, off_t orig_pos, MVPFile *m, MVPFile *m2);

/**
 * callback function for dct image hash use in mvptree structure.
 */
__declspec(dllexport)
float hammingdistance(DP *pntA, DP *pntB);

/** /brief aux function to query
 *  /param m - MVPFile state information
 *  /param query - DP of datapoint to query
 *  /param knearest - int capacity of results array.
 *  /param radius - float value of radius of values to consider
 *  /param results - DP array of points of result
 *  /param count - int* number of results found (out)
 *  /param level - int value to track recursion depth.
 *  /return MVPRetCode
**/
__declspec(dllexport)
MVPRetCode ph_query_mvptree(MVPFile *m, DP *query, int knearest, float radius,
			    DP **results, int *count, int level);

/**  /brief query mvptree function
 *   /param m - MVPFile file state info
 *   /param query - DP* item to query for
 *   /param knearest - int capacity of results array
 *   /param radius  - float radius to consider in query
 *   /param results - DP** list of pointers to results found
 *   /param count -  int number of results found (out)
 **/
__declspec(dllexport)
MVPRetCode ph_query_mvptree(MVPFile *m, DP *query, int knearest, float radius,
			    DP **results, int *count);

/** /brief save dp points to a file (aux func)
 *  /param m - MVPFile state information of file
 *  /param points - DP** list of points to add
 *  /param nbpoints - int length of points array
 *  /param saveall_flag - int  1 indicates used by save, 0 indicates used by add function
 *  /param level - int track recursion level
 *  /return FileIndex* - fileno and offset into file.
**/
__declspec(dllexport)
FileIndex* ph_save_mvptree(MVPFile *m, DP **points, int nbpoints, int saveall_flag, int level);

/** /brief save points to mvp file 
 *  /param m - MVPFile state info of file
 *  /param points - DP** list of points to add
 *  /param nbpoints - int number of points
 *  /return MVPRetCode - ret code 
**/
__declspec(dllexport)
MVPRetCode ph_save_mvptree(MVPFile *m, DP **points, int nbpoints);

/**  /brief add points to mvp file (aux function)
 *   /param m - MVPFile state information of file
 *   /param new_dp - datapoint to add
 *   /param level - int track recursion level
 *   /return MVPRetCode
 **/
__declspec(dllexport)
MVPRetCode ph_add_mvptree(MVPFile *m, DP *new_dp, int level);

/** /brief add a list of points to mvp file
    /param m - MVPFile state information of file.
    /param points - DP** list of points to add
    /param nbpoints - int number of points
    /return int - number of points added, neg for error
**/
__declspec(dllexport)
int ph_add_mvptree(MVPFile *m, DP **points, int nbpoints);

/** /brief textual hash for file
 *  /param filename - char* name of file
 *  /param nbpoints - int length of array of return value (out)
 *  /return TxtHashPoint* array of hash points with respective index into file.
 **/
__declspec(dllexport)
TxtHashPoint* ph_texthash(const char *filename, int *nbpoints);

/** /brief compare 2 text hashes
 *  /param hash1 -TxtHashPoint
 *  /param N1 - int length of hash1
 *  /param hash2 - TxtHashPoint
 *  /param N2 - int length of hash2
 *  /param nbmatches - int number of matches found (out)
 *  /return TxtMatch* - list of all matches
 **/
__declspec(dllexport)
TxtMatch* ph_compare_text_hashes(TxtHashPoint *hash1, int N1, TxtHashPoint *hash2, int N2, int *nbmatches);

#endif
