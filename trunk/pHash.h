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

#include <limits.h>
#include <math.h>
#include <dirent.h>
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



/* structure for a single image hash */
struct datapoint {
    char *id;
    ulong64 hash;
    int *path;
};

typedef struct datapoint DP;

/**   /brief struct for leaf node in mvptree
 *
 **/
struct leaf_node {
    int node_type;//0=leaf node, 1=internal_node
    DP *sv1, *sv2;
    DP **points;
    int *D1;
    int *D2;
    int Np;
};

typedef struct leaf_node Leaf;


union Node;

/**
 *   /brief struct for internal_node of mvptree
 **/
struct internal_node {
    int node_type; //0=leaf_node, 1=internal_node
    DP *sv1, *sv2;
    float *M1;
    float *M2;
    Node **child;
    int Nc;
};

typedef struct internal_node InternalNode;

//node of the mvptree
union Node {
    Leaf leaf;
    InternalNode internal;
};

typedef unsigned char uint8_t;

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


/**
 *  DO NOT CHANGE.  For use in experimenting with the size and shape of the mvp tree structure.
**/
#define M 2                      //number of brances in tree      (adjustable)
#define LENGTH_M1  M-1           //number of first level pivots   (leave alone)
#define LENGTH_M2  M*(LENGTH_M1) //number of 2nd level pivots     (leave alone)
#define FANOUT M*M               //fanout per node                (leave alone)
#define P 5                      //path length                    (adjustable)
#define K 25                     //number points, Np in leaf node (adjustable)
#define CAPACITY 500             //maximum number of files in tree(adjustable)

/* do not change */
static int BranchFactor = M;     //changed by save/read functions for mvptree
static int LengthM1 = LENGTH_M1;
static int LengthM2 = LENGTH_M2;
static int Fanout   = FANOUT;
static int PathLength = P;
static int LeafCapacity = K;


/* /brief alloc a single data point
 *  allocates path array, does nto set id or path
 */
DP* ph_malloc_datapoint();

/** /brief free a datapoint and its path
 *
 */
void ph_free_datapoint(DP *dp);

/**
 *    /brief alloc a leaf for an mvptree
 **/
Node* ph_malloc_leaf();

/**  /brief alloc an inner node for mvptree
 *
 **/
Node* ph_malloc_inner();

/** /brief free a node of an mvptree
 *
 **/
void ph_free_node(Node *node);

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

/** /brief creat an mvptree
 *  /param points - array of pointers to datapoint struct's
 *  /param nbpoints - int value number of datapoints
 *  /param level - int value to track recursion (user can ignore)
 *  /return Node pointer to top of created tree (NULL for error)
 **/
Node* ph_createMVPtree(DP **points, int nbpoints, int level=0);

/** /brief add a new datapoint to an existing tree
 *  /param tree - node pointer to top of an existing tree
 *  /param dp   - pointer to a new datapoint struct
 *  /param level - int value to track recursion (user can ignore)
 *  /return Node pointer to top of new tree (NULL for error)
 **/
Node* ph_addDPtoMVPtree(Node *tree,DP *dp,int level=0);

/** /brief create a list of datapoint's directly from a directory of image files
 *  /param dirname - path and name of directory containg all image file names
 *  /param capacity - int value for upper limit on number of hashes
 *  /param count - number of hashes created (out param)
 *  /return pointer to a list of DP pointers (NULL for error)
 **/
DP** ph_read_imagehashes(const char *dirname,int capacity, int &count);

/** /brief get all the filenames in specified directory
 *  /param dirname - string value for path and filename
 *  /param cap - int value for upper limit to number of files
 *  /param count - int value for number of file names returned
 *  /return array of pointers to string file names (NULL for error)
 **/
char** ph_readfilenames(const char *dirname,int cap,int &count);

/** /brief save a datapoint to file 
 *  auxiliary function for saving the mvp tree structure
 *  /param dp - pointer to a datapoint struct
 *  /param pfile - FILE pointer to an opened file
 *  /return 0 for success, -1 for error
 **/
int ph_saveDP(DP *dp, FILE *pfile);

/** /brief read a datapoint from file
 *  auxiliary function for reading mvp tree from file
 *  /param pfile - FILE pointer to opened file
 *  /return DP* - pointer value to a newly read datapoint, NULL for error
 **/
DP* ph_readDP(FILE *pfile);

/** /brief read mvp tree from file
 *  auxiliary function for reading from file
 *  /param pfile - FILE pointer to opened file
 *  /return Node ptr to top of tree, NULL for error
 **/
Node* ph_readMVPtree(FILE *pfile);

/** /brief read mvp tree from specified file
 *  /param filename - string value
 *  /return Node ptr to newly read tree, NULL on error
 **/ 
Node* ph_readMVPtree(const char *filename);

/** /brief save mvptree to file
 *  /param tree - Node ptr to tree
 *  /param pfile - FILE ptr to opened file
 *  /return int value, 0 for success, -1 for error 
 **/
int ph_saveMVPtree(Node *tree, FILE *pfile);

/** /brief save tree to specified file
 *  /param tree - Node ptr to tope of tree
 *  /param filename - string value 
 *  /return int value, 0 for success, -1 on error
 **/
int ph_saveMVPtree(Node *tree, const char *filename);

/** /brief print out a tree 
 *  for debugging purposes
 **/
void ph_printMVPtree(Node *tree);

/** /brief query mvp tree 
 *  /param tree - node ptr to top of tree 
 *  /param query - specific datapoint to query
 *  /param r - int value for radius of search query
 *  /param k - return k nearest neighbors
 *  /param path - int array of length P to track path distances from query to vantage points.(alloc before call)
 *  /param results - array of pointers to k datapoints (allocated before call)
 *  /param count - int value for number of results found (out param).
 *  /param level - int to track recursion depth (user can ignore).
 **/
int ph_queryMVPtree(Node *tree,DP *query,int r,int k,int *path,DP **results,int &count,int level=0);


#endif
