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

#ifndef _AUDIO_PHASH_H
#define _AUDIO_PHASH_H

#include <limits.h>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <algorithm>

extern "C" {
	#include "libavformat/avformat.h"
	#include "libavcodec/avcodec.h"
	#include "libswscale/swscale.h"
        #include "ph_fft.h"
}

/*  /brief count number of samples in file
 *  
 *  /param filename - path and file name of audio file
 *  /param sr - sample rate conversion
 *  /param channels - channels number conversion
 *  /return int count of number of sampels, negative for error
*/
int ph_count_samples(const char *filename, int sr,int channels);




/* /brief read audio 
 *
 * /param filename - path and name of audio file to read
 * /param sr - sample rate conversion
 * /param channels - number channel conversion
 * /param N - (out) param for buf length
 * /param float* - float pointer to start of buffer, NULL if error
 */
float* ph_readaudio(const char *filename, int sr, int channels, int &N);



/* /brief audio hash calculation
 * purpose: hash calculation for each frame in the buffer.
 *          Each value is computed from successive overlapping frames of the input buffer. 
 *          The value is based on the bark scale values of the frame fft spectrum. The value
 *          computed from temporal and spectral differences on the bark scale.
 * 
 * /param buf - pointer to start of buffer
 * /param N   - length of buffer
 * /param sr  - sample rate on which to base the audiohash
 * /param nb_frames - (out) number of frames in audio buf and length of audiohash buffer returned
 * /return uint32 pointer to audio hash, NULL for error
*/
uint32_t* ph_audiohash(float *buf, int N, int sr, int &nb_frames);



/* /brief bit count set bits in 32bit variable
 * /param n 
 * /return int number of bits set to 1, negative if error
 */
int ph_bitcount(uint32_t n);


/* /brief compare 2 hash blocks 
 * /param ptr_blockA - pointer to the first block
 * /param ptr_blockB - pointer to the second block
 * /param block_size - length of both blocks to compare
 * /return double bit error rate (ber) from comparing two blocks, neg for error
 */
double ph_compare_blocks(const uint32_t *ptr_blockA,const uint32_t *ptr_blockB, const int block_size);


/* /brief distance function between two hashes
 *
 * /param hash_a - first hash
 * /param Na     - length of first hash
 * /param hash_b - second hash
 * /param Nb     - length of second hash
 * /param threshold - threshold value to compare successive blocks, 0.25, 0.30, 0.35
 * /param block_size - length of block_size, 256
 * /param Nc     - (out) length of confidence score vector
 * /return double - ptr to confidence score vector
 */
double* ph_audio_distance_ber(uint32_t *hash_a , const int Na, uint32_t *hash_b, const int Nb, const float threshold, const int block_size, int &Nc);


#endif
