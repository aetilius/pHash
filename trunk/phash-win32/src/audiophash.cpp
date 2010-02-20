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

#include <limits.h>
#include <math.h>
#include "audiophash.h"

extern "C" {
	#include "libavformat/avformat.h"
	#include "libavcodec/avcodec.h"
	#include "libswscale/swscale.h"
    #include "phcomplex.h"
    #include "ph_fft.h"
}

__declspec(dllexport)
int ph_count_samples(const char *filename, int sr,int channels){

	av_register_all();
	avcodec_register_all();

	AVFormatContext *pFormatCtx;
	
	// Open file
	if(av_open_input_file(&pFormatCtx, filename, NULL, 0, NULL)!=0)
	  return -1 ; // Couldn't open file
	 
	// Retrieve stream information
	if(av_find_stream_info(pFormatCtx)<0)
	  return -1; // Couldn't find stream information
	
	//dump_format(pFormatCtx,0,NULL,0);//debugging function to print infomation about format

	unsigned int i;
	AVCodecContext *pCodecCtx;

	// Find the video stream
	int audioStream=-1;

	for(i=0; i<pFormatCtx->nb_streams; i++)
	{
	     if(pFormatCtx->streams[i]->codec->codec_type==CODEC_TYPE_AUDIO) 
	     {
		    audioStream=i;
		    break;
	     }
	}
	if(audioStream==-1){
	    printf("no audio stream\n");
	     return -1; //no video stream
	}
	
	// Get a pointer to the codec context for the audio stream
	pCodecCtx=pFormatCtx->streams[audioStream]->codec;

	int src_sr = pCodecCtx->sample_rate;
        int src_channels = pCodecCtx->channels;

	AVCodec *pCodec;

	// Find the decoder
	pCodec=avcodec_find_decoder(pCodecCtx->codec_id);
	if(pCodec==NULL) 
	{
	  	return -1 ; // Codec not found
	}
	// Open codec
	if(avcodec_open(pCodecCtx, pCodec)<0)
	  return -1; // Could not open codec

	uint8_t *in_buf = (uint8_t *)malloc(AVCODEC_MAX_AUDIO_FRAME_SIZE + FF_INPUT_BUFFER_PADDING_SIZE);
        uint8_t *out_buf = (uint8_t *)malloc(AVCODEC_MAX_AUDIO_FRAME_SIZE);

        int in_buf_used, numbytesread, buf_size = AVCODEC_MAX_AUDIO_FRAME_SIZE;


	ReSampleContext *rs_ctxt = audio_resample_init(channels, src_channels, sr, src_sr);

	int count = 0;
	AVPacket packet;
	while(av_read_frame(pFormatCtx, &packet)>=0) 
	{
	    in_buf_used = buf_size;
	    numbytesread = avcodec_decode_audio2(pCodecCtx,(int16_t*)in_buf,&in_buf_used,packet.data,packet.size);  
	    if (numbytesread <= 0){
		fprintf(stderr,"error reading from audiostream\n");
		continue;
	    }
	    count += (int)(audio_resample(rs_ctxt,(short*)out_buf ,(short*)in_buf,in_buf_used)/sizeof(int16_t));
	}
	free(in_buf);
	free(out_buf);
	audio_resample_close(rs_ctxt);
	avcodec_close(pCodecCtx);
	av_close_input_file(pFormatCtx);
	
	return count;
}

__declspec(dllexport)
float* ph_readaudio(const char *filename, int sr, int channels, float *sigbuf, int &buflen)
{
	int cap = 0;

	av_register_all();

	AVFormatContext *pFormatCtx;
	
	// Open file
	if(av_open_input_file(&pFormatCtx, filename, NULL, 0, NULL)!=0){
	    buflen=0;
	  return NULL ; // Couldn't open file
	}
	 
	// Retrieve stream information
	if(av_find_stream_info(pFormatCtx)<0){
	    buflen=0;
	    return NULL; // Couldn't find stream information
	}
	//dump_format(pFormatCtx,0,NULL,0);//debugging function to print infomation about format

	unsigned int i;
	AVCodecContext *pCodecCtx;

	// Find the video stream
	int audioStream=-1;

	for(i=0; i<pFormatCtx->nb_streams; i++)
	{
	     if(pFormatCtx->streams[i]->codec->codec_type==CODEC_TYPE_AUDIO) 
	     {
		    audioStream=i;
		    break;
	     }
	}
	if(audioStream==-1){
	     buflen = 0;
	     return NULL; //no video stream
	}
    
	// Get a pointer to the codec context for the audio stream
	pCodecCtx=pFormatCtx->streams[audioStream]->codec;

    int duration = (int)pFormatCtx->streams[audioStream]->duration;
	int src_sr = pCodecCtx->sample_rate;
    int src_channels = pCodecCtx->channels;
    
    float *buf; 
	if ((!sigbuf)||(buflen <= duration)){
        //buf = (float*)av_malloc(duration*sizeof(float)); /* alloc buffer */
        //buflen = (int)duration;
        buflen = 0; /* buffer too small, bail out */
        return NULL;
	}
	else {
        buf = sigbuf;   /* use buffer handed in param */     
	}
	AVCodec *pCodec;

	// Find the decoder
	pCodec=avcodec_find_decoder(pCodecCtx->codec_id);
	if(pCodec==NULL) 
	{
	    buflen=0;
	  	return NULL ; // Codec not found
	}
	// Open codec
	if(avcodec_open(pCodecCtx, pCodec)<0){
	    buflen=0;
	    return NULL; // Could not open codec
	}

	uint8_t *in_buf = (uint8_t*)av_malloc(AVCODEC_MAX_AUDIO_FRAME_SIZE + FF_INPUT_BUFFER_PADDING_SIZE);
    int in_buf_used, numbytesread, buf_size = AVCODEC_MAX_AUDIO_FRAME_SIZE;
	uint8_t *out_buf = (uint8_t*)av_malloc(AVCODEC_MAX_AUDIO_FRAME_SIZE);
    int16_t *out_buf16 = (int16_t*)out_buf;

	ReSampleContext *rs_ctx = audio_resample_init(channels, src_channels, sr, src_sr);
    int index = 0;
	AVPacket packet;
	while(av_read_frame(pFormatCtx, &packet)>=0) 
	{
            
	    in_buf_used = buf_size;
	    numbytesread = avcodec_decode_audio2(pCodecCtx,(int16_t*)in_buf,&in_buf_used,packet.data,packet.size);  
	    if (numbytesread <= 0){
		    continue;
	    }
	    int cnt = audio_resample(rs_ctx,(short*)out_buf,(short*)in_buf,(int)(in_buf_used/sizeof(int16_t)));
		if (index + cnt > buflen){ /*exceeds capacity of buffer, bail out */ 
            buf = NULL;
            buflen = 0;
            goto readaudio_cleanup;
		}	
		for (int i=0;i<cnt;i++){
           buf[index+i] = (float)out_buf16[i]/(float)SHRT_MAX;
		}

		cap   += cnt;
		index += cnt;
	}

readaudio_cleanup:
	av_free(in_buf);
	av_free(out_buf);
	audio_resample_close(rs_ctx);
	avcodec_close(pCodecCtx);
	av_close_input_file(pFormatCtx);
    buflen = cap;
	
	return buf;
} 

__declspec(dllexport)
uint32_t* ph_audiohash(const float *buf, const int nbbuf, uint32_t *hashbuf, const int nbcap, const int sr,int &nbframes){

   const int frame_length = 4096;//2^12
   int nfft = frame_length;
   const int nfft_half = 2048;
   int start = 0;
   int end = start + frame_length - 1;
   int overlap = (int)(31*frame_length/32);
   int advance = frame_length - overlap;
   int index = 0;
   nbframes = (int)(floor(float(nbbuf/advance)) - floor(float(frame_length/advance)) + 1);
   double window[frame_length];
   for (int i = 0;i<frame_length;i++){
       //hamming window
       window[i] = 0.54 - 0.46*cos(2*M_PI*i/(frame_length-1));
   }
   
   double frame[frame_length];
   Complexd *pF = (Complexd*)malloc(sizeof(Complexd)*nfft);

   double magnF[nfft_half];
   double maxF = 0.0;
   double maxB = 0.0;
  
   double minfreq = 300;
   double maxfreq = 3000;
   double temp = double(minfreq/600.0);
   double minbark = 6*log(temp + sqrt(temp*temp + 1.0));  
   temp = double(maxfreq/600.0);
   double maxbark = 6*log(temp + sqrt(temp*temp + 1.0));
   double nyqbark = maxbark - minbark;
   const int nfilts = 33;
   double stepbarks = nyqbark/(nfilts - 1);
   int nb_barks = (int)(floor(float(nfft_half/2 + 1)));
   double barkwidth = 1.06;    

   double *freqs = new double[nb_barks];
   double *binbarks = new double[nb_barks];
   double curr_bark[nfilts];
   double prev_bark[nfilts];
   for (int i=0;i< nfilts;i++){
       prev_bark[i] = 0.0;
   }
   uint32_t *hash = NULL;
   if (!hashbuf || (nbframes > nbcap)){
	   //hash = (uint32_t*)malloc(nbframes*sizeof(uint32_t)); /*must allocate new buffer */
       free(pF);
       delete [] freqs;
       delete [] binbarks;
       nbframes = 0; /*hash buf not big enough, bail out */
       return NULL;
   }
   else {
       hash = hashbuf; /*use buf provided in parameters */
   }

   double lof,hif;

   for (int i=0; i < nb_barks;i++){
	   temp = i*sr/nfft_half/600.0;
       binbarks[i] = 6*log(temp + sqrt(temp*temp + 1.0));
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
           double m = ph_min(lof,hif);
           m = ph_min(0.0,m);
           m = pow(10,m);
           wts[i][j] = m;
       }
   }

   while (end <= nbbuf){
       maxF = 0.0;
       maxB = 0.0;
       for (int i = 0;i<frame_length;i++){
		   frame[i] = window[i]*buf[start+i];
       }
       fft(frame, frame_length, pF);

       for (int i=0; i < nfft_half;i++){
		   magnF[i] = sqrt(pF[i].re*pF[i].re + pF[i].im*pF[i].im);
		   if (magnF[i] > maxF)
			   maxF = magnF[i];
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

   free(pF);
   delete [] freqs;
   delete [] binbarks;
   for (int i=0;i<nfilts;i++){
       delete [] wts[i];
   }
   delete [] wts;

   return hash;
}
__declspec(dllexport)
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
__declspec(dllexport)
double ph_compare_blocks(const uint32_t *ptr_blockA,const uint32_t *ptr_blockB, const int block_size){
    double result = 0.0f;
    for (int i=0;i<block_size;i++){
		uint32_t xordhash = ptr_blockA[i]^ptr_blockB[i];
        result += ph_bitcount(xordhash);
    }
    result = result/(32*block_size);
    return result;
}
__declspec(dllexport)
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

    double *pC=NULL;
	if ((N1 < block_size) || (N2 < block_size)){
        Nc = 0;
        return NULL;
	} 
    pC = new double[Nc];
	if (pC == NULL){
        Nc = 0;
        return NULL;
	}
    int M = (int)((float)N1/(float)block_size);
    double *dist = new double[M];
	if (!dist){
         delete pC;
         Nc = 0;
         return NULL;
	}
    for (int i=0; i < Nc;i++){
        uint32_t *pha = ptrA;
        uint32_t *phb = ptrB + i;
        int k=0;
		while (k < M){	
			dist[k++] = ph_compare_blocks(pha,phb,block_size);
            pha += block_size;
            phb += block_size;
		};
		double sum_above = 0.0f;
		double sum_below = 0.0f;
		int nb_above = 0;
		int nb_below = 0;
		for (int n = 0; n < M; n++){
			if (dist[n] <= threshold){
				sum_below += 1-dist[n];
				nb_below++;
			} else {
				sum_above += 1-dist[n];
				nb_above++;
			}
		}
		double above_factor = sum_above/(double)M;
		double below_factor = sum_below/(double)M;
		pC[i] = 0.5*(1 + below_factor - above_factor);

	}
    delete [] dist;

    return pC;
}
