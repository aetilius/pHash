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

float* ph_readaudio(const char *filename, int sr, int channels, int &N)
{
	N = 500000;
	int cap = 0;
	float *buf = (float*)malloc(N*sizeof(float));


	av_register_all();

	AVFormatContext *pFormatCtx;
	
	// Open file
	if(av_open_input_file(&pFormatCtx, filename, NULL, 0, NULL)!=0){
	    N=0;
	  return NULL ; // Couldn't open file
	}
	 
	// Retrieve stream information
	if(av_find_stream_info(pFormatCtx)<0){
	    N=0;
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
	     N = 0;
	     return NULL; //no video stream
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
	        N=0;
	  	return NULL ; // Codec not found
	}
	// Open codec
	if(avcodec_open(pCodecCtx, pCodec)<0){
	    N=0;
	    return NULL; // Could not open codec
	}

	uint8_t *in_buf = (uint8_t*)malloc(AVCODEC_MAX_AUDIO_FRAME_SIZE + FF_INPUT_BUFFER_PADDING_SIZE);
        int in_buf_used, numbytesread, buf_size = AVCODEC_MAX_AUDIO_FRAME_SIZE;
	uint8_t *out_buf = (uint8_t*)malloc(AVCODEC_MAX_AUDIO_FRAME_SIZE);
        int16_t *out_buf16 = (int16_t*)out_buf;

	ReSampleContext *rs_ctx = audio_resample_init(channels, src_channels, sr, src_sr);
        int index = 0;
	AVPacket packet;
	while(av_read_frame(pFormatCtx, &packet)>=0) 
	{
            
	    in_buf_used = buf_size;
	    numbytesread = avcodec_decode_audio2(pCodecCtx,(int16_t*)in_buf,&in_buf_used,packet.data,packet.size);  
	    if (numbytesread <= 0){
		fprintf(stderr,"error reading from audiostream\n");
		continue;
	    }
	    int cnt = audio_resample(rs_ctx,(short*)out_buf,(short*)in_buf,(int)(in_buf_used/sizeof(int16_t)));
			
	    if (cap + cnt > N){
		float *tmpbuf = (float*)realloc(buf,(N + 2*cnt)*sizeof(float));
		if (!tmpbuf){
		    N=0;
		    return NULL;
		}
		buf = tmpbuf;
                N += 2*cnt;
	    }   
	   for (int i=0;i<cnt;i++){
	       buf[index+i] = ((float)out_buf16[i]/(float)SHRT_MAX);

	   }

           cap   += cnt;
           index += cnt;
	}

	free(in_buf);
	free(out_buf);
	audio_resample_close(rs_ctx);
	avcodec_close(pCodecCtx);
	av_close_input_file(pFormatCtx);
        N = cap;
	
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
   fftw_complex *pF;
   fftw_plan p;

   pF = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nfft);
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

   p = fftw_plan_dft_r2c_1d(frame_length,frame,pF,FFTW_ESTIMATE);

   while (end <= N){
       maxF = 0.0;
       maxB = 0.0;
       for (int i = 0;i<frame_length;i++){
	   frame[i] = window[i]*buf[start+i];
       }
       fftw_execute(p);
       for (int i=0; i < nfft_half;i++){
	   magnF[i] = sqrt(pF[i][0]*pF[i][0] +  pF[i][1]*pF[i][1] );
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
   


   fftw_destroy_plan(p);
   fftw_free(pF);
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
