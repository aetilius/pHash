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

#include "audiophash.h"
#include <sndfile.h>
#include <samplerate.h>

int ph_count_samples(const char *filename, int sr,int channels){

    SF_INFO sf_info;
    sf_info.format=0;
    SNDFILE *sndfile = sf_open(filename, SFM_READ, &sf_info);
    if (sndfile == NULL){
	return NULL;
    }
    int count = sf_info.frames;
    sf_close(sndfile);
    return count;
}

float* ph_readaudio(const char *filename, int sr, int channels, float *sigbuf, int &buflen, const float nbsecs)
{
    SF_INFO sf_info;
    sf_info.format=0;
    SNDFILE *sndfile = sf_open(filename, SFM_READ, &sf_info);
    if (sndfile == NULL){
	return NULL;
    }    

    sf_command(sndfile, SFC_SET_NORM_FLOAT, NULL, SF_TRUE);

    //allocate input buffer for signal
    int src_frames = (nbsecs <= 0) ? sf_info.frames : (int)(nbsecs*sf_info.samplerate);
    src_frames = (sf_info.frames < src_frames) ? sf_info.frames : src_frames;
    float *inbuf = (float*)malloc(src_frames*sf_info.channels*sizeof(float));
    if (!inbuf){
	sf_close(sndfile);
	return NULL;
    }

    sf_count_t cnt_frames = sf_readf_float(sndfile, inbuf, src_frames);
    
    double sr_ratio = (double)sr/(double)sf_info.samplerate;
    if (src_is_valid_ratio(sr_ratio) == 0){
	sf_close(sndfile);
	free(inbuf);
	return NULL;
    }

    int inbuflen = cnt_frames*sf_info.channels;   
    int outbuflen = (int)(sr_ratio*cnt_frames*sf_info.channels);
    float *outbuf = (float*)malloc(outbuflen*sizeof(float));
    if (!outbuf){
	sf_close(sndfile);
	free(inbuf);
	return NULL;
    }

    int error;
    SRC_STATE *src_state = src_new(SRC_LINEAR, sf_info.channels, &error);
    if (!src_state){
	sf_close(sndfile);
	free(inbuf);
	free(outbuf);
        return NULL;
    }

    SRC_DATA src_data;
    src_data.data_in = inbuf;
    src_data.data_out = outbuf;
    src_data.input_frames = cnt_frames;
    src_data.output_frames = outbuflen/sf_info.channels;
    src_data.end_of_input = SF_TRUE;
    src_data.src_ratio = sr_ratio;

    if (error = src_process(src_state, &src_data)){
        sf_close(sndfile);
	free(inbuf);
	free(outbuf);
	src_delete(src_state);
	return NULL;
    }

    float *buf = NULL;
    if ((buflen > src_data.output_frames*channels) || (sigbuf == NULL)){
	//alloc new buffer
	buf = (float*)malloc(src_data.output_frames*channels*sizeof(float));
    } else {
	//use buffer from param
	buf = sigbuf;
    }
    buflen = src_data.output_frames*channels;

    if (!buf){
	sf_close(sndfile);
	free(inbuf);
	free(outbuf);
	src_delete(src_state);
	return NULL;
    }
    buflen = src_data.output_frames;

    if (channels == 1){
	//average across all channels
	for (int i=0;i<src_data.output_frames*sf_info.channels;i+=sf_info.channels){
	    buf[i/sf_info.channels] = 0;
	    for (int j=0;j<sf_info.channels;j++){
		buf[i/sf_info.channels] += outbuf[i+j];
	    }
	    buf[i/sf_info.channels] /= sf_info.channels;
	}
    } else if (channels <= sf_info.channels){
	//just grab first nb channels
	for (int i=0;i<src_data.output_frames*sf_info.channels;i+=sf_info.channels){
	    for (int j=0;j<channels;j++){
		buf[i/sf_info.channels+j] = outbuf[i+j];
	    }
	}
    } else {
	free(buf);
	buf=NULL;
    }

    src_delete(src_state);
    sf_close(sndfile);
    free(inbuf);
    free(outbuf);
    return buf;
} 

uint32_t* ph_audiohash(float *buf, int N, int sr, int &nb_frames){

   int frame_length = 4096;//2^12
   int nfft = frame_length;
   int nfft_half = 2048;
   int start = 0;
   int end = start + frame_length - 1;
   int overlap = (int)(31*frame_length/32);
   int advance = frame_length - overlap;
   int index = 0;
   nb_frames = (int)(floor(N/advance) - floor(frame_length/advance) + 1);
   double window[frame_length];
   for (int i = 0;i<frame_length;i++){
       //hamming window
       window[i] = 0.54 - 0.46*cos(2*M_PI*i/(frame_length-1));
   }
   
   double frame[frame_length];
   //fftw_complex *pF;
   //fftw_plan p;
   //pF = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nfft);
   complex double *pF = (complex double*)malloc(sizeof(complex double)*nfft);

   double magnF[nfft_half];
   double maxF = 0.0;
   double maxB = 0.0;
  
   double minfreq = 300;
   double maxfreq = 3000;
   double minbark = 6*asinh(minfreq/600.0);
   double maxbark = 6*asinh(maxfreq/600.0);
   double nyqbark = maxbark - minbark;
   int nfilts = 33;
   double stepbarks = nyqbark/(nfilts - 1);
   int nb_barks = (int)(floor(nfft_half/2 + 1));
   double barkwidth = 1.06;    

   double freqs[nb_barks];
   double binbarks[nb_barks];
   double curr_bark[nfilts];
   double prev_bark[nfilts];
   for (int i=0;i< nfilts;i++){
       prev_bark[i] = 0.0;
   }
   uint32_t *hash = (uint32_t*)malloc(nb_frames*sizeof(uint32_t));
   double lof,hif;

   for (int i=0; i < nb_barks;i++){
       binbarks[i] = 6*asinh(i*sr/nfft_half/600.0);
       freqs[i] = i*sr/nfft_half;
   }
   double **wts = new double*[nfilts];
   for (int i=0;i<nfilts;i++){
      wts[i] = new double[nfft_half];
   }
   for (int i=0;i<nfilts;i++){
       for (int j=0;j<nfft_half;j++){
	   wts[i][j] = 0.0;
       }
   }
  
   //calculate wts for each filter
   for (int i=0;i<nfilts;i++){
       double f_bark_mid = minbark + i*stepbarks;
       for (int j=0;j<nb_barks;j++){
	   double barkdiff = binbarks[j] - f_bark_mid;
           lof = -2.5*(barkdiff/barkwidth - 0.5);
           hif = barkdiff/barkwidth + 0.5;
           double m = std::min(lof,hif);
           m = std::min(0.0,m);
           m = pow(10,m);
           wts[i][j] = m;
       }
   }

   //p = fftw_plan_dft_r2c_1d(frame_length,frame,pF,FFTW_ESTIMATE);

   while (end <= N){
       maxF = 0.0;
       maxB = 0.0;
       for (int i = 0;i<frame_length;i++){
	   frame[i] = window[i]*buf[start+i];
       }
       //fftw_execute(p);
       if (fft(frame, frame_length, pF) < 0){
	   return NULL;
       }
       for (int i=0; i < nfft_half;i++){
	   //magnF[i] = sqrt(pF[i][0]*pF[i][0] +  pF[i][1]*pF[i][1] );
           magnF[i] = cabs(pF[i]);
	   if (magnF[i] > maxF){
	       maxF = magnF[i];
	   }
       }

       for (int i=0;i<nfilts;i++){
	   curr_bark[i] = 0;
	   for (int j=0;j < nfft_half;j++){
	       curr_bark[i] += wts[i][j]*magnF[j];
	   }
           if (curr_bark[i] > maxB)
	       maxB = curr_bark[i];
       }

       uint32_t curr_hash = 0x00000000u;
       for (int m=0;m<nfilts-1;m++){
	   double H = curr_bark[m] - curr_bark[m+1] - (prev_bark[m] - prev_bark[m+1]);
	   curr_hash = curr_hash << 1;
	   if (H > 0)
	       curr_hash |= 0x00000001;
       }


       hash[index] = curr_hash;
       for (int i=0;i<nfilts;i++){
	   prev_bark[i] = curr_bark[i];
       }
       index += 1;
       start += advance;
       end   += advance;
   }
   


   //fftw_destroy_plan(p);
   //fftw_free(pF);
   free(pF);
   for (int i=0;i<nfilts;i++){
       delete [] wts[i];
   }
   delete [] wts;
   return hash;
}


int ph_bitcount(uint32_t n){
    
    //parallel bit count
    #define MASK_01010101 (((uint32_t)(-1))/3)
    #define MASK_00110011 (((uint32_t)(-1))/5)
    #define MASK_00001111 (((uint32_t)(-1))/17)

    n = (n & MASK_01010101) + ((n >> 1) & MASK_01010101) ;
    n = (n & MASK_00110011) + ((n >> 2) & MASK_00110011) ;
    n = (n & MASK_00001111) + ((n >> 4) & MASK_00001111) ;
    return n % 255;

}

double ph_compare_blocks(const uint32_t *ptr_blockA,const uint32_t *ptr_blockB, const int block_size){
    double result = 0;
    for (int i=0;i<block_size;i++){
	uint32_t xordhash = ptr_blockA[i]^ptr_blockB[i];
        result += ph_bitcount(xordhash);
    }
    result = result/(32*block_size);
    return result;
}
double* ph_audio_distance_ber(uint32_t *hash_a , const int Na, uint32_t *hash_b, const int Nb, const float threshold, const int block_size, int &Nc){

    uint32_t *ptrA, *ptrB;
    int N1, N2;
    if (Na <= Nb){
	ptrA = hash_a;
	ptrB = hash_b;
	Nc = Nb - Na + 1;
	N1 = Na;
        N2 = Nb;
    } else {
	ptrB = hash_a;
	ptrA = hash_b;
	Nc = Na - Nb + 1;
	N1 = Nb;
	N2 = Na;
    }

    double *pC = new double[Nc];
    if (!pC)
	return NULL;
    int k,M,nb_above, nb_below, hash1_index,hash2_index;
    double sum_above, sum_below,above_factor, below_factor;

    uint32_t *pha,*phb;
    double *dist = NULL;

    for (int i=0; i < Nc;i++){

	M = (int)floor(std::min(N1,N2-i)/block_size);

        pha = ptrA;
        phb = ptrB + i;

	double *tmp_dist = (double*)realloc(dist, M*sizeof(double));
        if (!tmp_dist){
	    return NULL;
        }
        dist = tmp_dist;
	dist[0] = ph_compare_blocks(pha,phb,block_size);

	k = 1;

	pha += block_size;
	phb += block_size;

	hash1_index = block_size;
	hash2_index = i + block_size;

	while ((hash1_index < N1 - block_size)  && (hash2_index < N2 - block_size)){
	    dist[k++] = ph_compare_blocks(pha,phb,block_size);
	    hash1_index += block_size;
	    hash2_index += block_size;
	    pha += block_size;
	    phb += block_size;
	}
        sum_above = 0;
	sum_below = 0;
	nb_above = 0;
	nb_below = 0;
	for (int n = 0; n < M; n++){

	    if (dist[n] <= threshold){
		sum_below += 1-dist[n];
		nb_below++;
	    } else {
		sum_above += 1-dist[n];
		nb_above++;
	    }
	}
	above_factor = sum_above/M;
	below_factor = sum_below/M;
	pC[i] = 0.5*(1 + below_factor - above_factor);
    }

    free(dist);
    return pC;
}
#ifdef HAVE_PTHREAD

void *ph_audio_thread(void *p)
{
        slice *s = (slice *)p;
        for(int i = 0; i < s->n; ++i)
        {
                DP *dp = (DP *)s->hash_p[i];
                int N, count;
		pair<int,int> *p = (pair<int,int> *)s->hash_params;
                float *buf = ph_readaudio(dp->id, p->first, p->second, NULL, N);
                uint32_t *hash = ph_audiohash(buf, N, p->first, count);
                free(buf);
                buf = NULL;
                dp->hash = hash;
                dp->hash_length = count;
        }
}

DP** ph_audio_hashes(char *files[], int count, int sr, int channels, int threads)
{
        if(!files || count == 0)
	        return NULL;

        int num_threads;
        if(threads > count)
        {
                num_threads = count;
        }
	else if(threads > 0)
        {
                num_threads = threads;
        }
	else
	{
                num_threads = ph_num_threads();
        }

	DP **hashes = (DP**)malloc(count*sizeof(DP*));

        for(int i = 0; i < count; ++i)
        {
                hashes[i] = (DP *)malloc(sizeof(DP));
                hashes[i]->id = strdup(files[i]);
        }

	pthread_t thds[num_threads];
	
        int rem = count % num_threads;
        int start = 0;
        int off = 0;
        slice *s = new slice[num_threads];
        for(int n = 0; n < num_threads; ++n)
        {
                off = (int)floor((count/(float)num_threads) + (rem>0?num_threads-(count % num_threads):0));

                s[n].hash_p = &hashes[start];
                s[n].n = off;
		s[n].hash_params = new pair<int,int>(sr,channels);
                start = off;
                --rem;
                pthread_create(&thds[n], NULL, ph_audio_thread, &s[n]);
        }
	for(int i = 0; i < num_threads; ++i)
        {
                pthread_join(thds[i], NULL);
		delete (pair<int,int>*)s[i].hash_params;
        }
        delete[] s;

	return hashes;

}
#endif
