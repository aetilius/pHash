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
    D Grant Starkweather - dstarkweather@phash.org

*/


#include "cimgffmpeg.h"

void vfinfo_close(VFInfo  *vfinfo){
    if (vfinfo->pFormatCtx != NULL){
	avcodec_close(vfinfo->pCodecCtx);
	vfinfo->pCodecCtx = NULL;
	avformat_close_input(&vfinfo->pFormatCtx);
	vfinfo->pFormatCtx = NULL;
	vfinfo->width = -1;
	vfinfo->height = -1;
    }
}

int ReadFrames(VFInfo *st_info, CImgList<uint8_t> *pFrameList, unsigned int low_index, unsigned int hi_index)
{
        //target pixel format
	PixelFormat ffmpeg_pixfmt;
	if (st_info->pixelformat == 0)
	    ffmpeg_pixfmt = PIX_FMT_GRAY8;
	else 
	    ffmpeg_pixfmt = PIX_FMT_RGB24;

	st_info->next_index = low_index;

	if (st_info->pFormatCtx == NULL){
	    st_info->current_index= 0;

        av_log_set_level(AV_LOG_QUIET);
	    av_register_all();
	
	    // Open video file
	    if(avformat_open_input(&st_info->pFormatCtx, st_info->filename, NULL, NULL)!=0)
		return -1 ; // Couldn't open file
	 
	    // Retrieve stream information
	    if(avformat_find_stream_info(st_info->pFormatCtx,NULL)<0)
		return -1; // Couldn't find stream information
	
	    //dump_format(pFormatCtx,0,NULL,0);//debugging function to print infomation about format
	
	    unsigned int i;
	    // Find the video stream
	    for(i=0; i<st_info->pFormatCtx->nb_streams; i++)
	    {
		if(st_info->pFormatCtx->streams[i]->codec->codec_type==AVMEDIA_TYPE_VIDEO) 
	        {
		    st_info->videoStream=i;
		    break;
		}
	    }
	    if(st_info->videoStream==-1)
		return -1; //no video stream
	
	
	    // Get a pointer to the codec context for the video stream
	    st_info->pCodecCtx = st_info->pFormatCtx->streams[st_info->videoStream]->codec;
	    if (st_info->pCodecCtx == NULL){
		return -1;
	    }

	    // Find the decoder
	    st_info->pCodec = avcodec_find_decoder(st_info->pCodecCtx->codec_id);
	    if(st_info->pCodec==NULL) 
	    {
	  	return -1 ; // Codec not found
	    }
	    // Open codec
	    if(avcodec_open2(st_info->pCodecCtx, st_info->pCodec,NULL)<0)
		return -1; // Could not open codec

	    st_info->height = (st_info->height<=0) ? st_info->pCodecCtx->height : st_info->height;
	    st_info->width  = (st_info->width<= 0) ? st_info->pCodecCtx->width : st_info->width;
	}

        AVFrame *pFrame;

	// Allocate video frame
	pFrame=avcodec_alloc_frame();
	if (pFrame==NULL)
	    return -1;

	// Allocate an AVFrame structure
	AVFrame *pConvertedFrame = avcodec_alloc_frame();
	if(pConvertedFrame==NULL)
	  return -1;
		
	uint8_t *buffer;
	int numBytes;
	// Determine required buffer size and allocate buffer
	numBytes=avpicture_get_size(ffmpeg_pixfmt, st_info->width,st_info->height);
	buffer=(uint8_t *)av_malloc(numBytes*sizeof(uint8_t));
	if (buffer == NULL)
	    return -1;

	avpicture_fill((AVPicture *)pConvertedFrame,buffer,ffmpeg_pixfmt,st_info->width,st_info->height);
		
	int frameFinished;
	int size = 0;
	

        int channels = ffmpeg_pixfmt == PIX_FMT_GRAY8 ? 1 : 3;

	AVPacket packet;
	int result = 1;
	CImg<uint8_t> next_image;
	SwsContext *c = sws_getContext(st_info->pCodecCtx->width, st_info->pCodecCtx->height, st_info->pCodecCtx->pix_fmt, st_info->width, st_info->height, ffmpeg_pixfmt , SWS_BICUBIC, NULL, NULL, NULL);
	while ((result>=0)&&(size<st_info->nb_retrieval)&&(st_info->current_index<=hi_index)){  
	  result =  av_read_frame(st_info->pFormatCtx, &packet);
          if (result < 0)
	      break;
    	  if(packet.stream_index==st_info->videoStream) {

		AVPacket avpkt; 
		av_init_packet(&avpkt); 
		avpkt.data = packet.data; 
		avpkt.size = packet.size; 
		// 
		// HACK for CorePNG to decode as normal PNG by default 
		// same method used by ffmpeg 
		avpkt.flags = AV_PKT_FLAG_KEY; 
	      
	 	avcodec_decode_video2(st_info->pCodecCtx, pFrame, &frameFinished,&avpkt);

	      // avcodec_decode_video(st_info->pCodecCtx, pFrame, &frameFinished,packet.data, packet.size);

	      if(frameFinished) {
		  if (st_info->current_index == st_info->next_index){
		      st_info->next_index += st_info->step;
		      sws_scale(c, pFrame->data, pFrame->linesize, 0, st_info->pCodecCtx->height, pConvertedFrame->data, pConvertedFrame->linesize);
			
		  next_image.assign(*pConvertedFrame->data, channels,st_info->width,st_info->height,1,true);
		  next_image.permute_axes("yzcx");
		  pFrameList->push_back(next_image);
		  size++;
		  }    
		  st_info->current_index++;
	      }
	      av_free_packet(&packet);
	  }
	}


	if (result < 0){
	    avcodec_close(st_info->pCodecCtx);
	    avformat_close_input(&st_info->pFormatCtx);
	    st_info->pFormatCtx = NULL;
	    st_info->pCodecCtx = NULL;
	    st_info->width = -1;
	    st_info->height = -1;
	}

	av_free(buffer);
	buffer = NULL;
	av_free(pConvertedFrame);
	pConvertedFrame = NULL;
	av_free(pFrame);
	pFrame = NULL;
	sws_freeContext(c);
	c = NULL;

	return size; 
}


int NextFrames(VFInfo *st_info, CImgList<uint8_t> *pFrameList)
{
        PixelFormat ffmpeg_pixfmt;
	if (st_info->pixelformat == 0)
	    ffmpeg_pixfmt = PIX_FMT_GRAY8;
        else 
	    ffmpeg_pixfmt = PIX_FMT_RGB24;

	if (st_info->pFormatCtx == NULL)
	{
	        st_info->current_index = 0;
		st_info->next_index = 0;
		av_register_all();
		st_info->videoStream = -1;
		//st_info->pFormatCtx = (AVFormatContext*)malloc(sizeof(AVFormatContext));
		//st_info->pCodecCtx = (AVCodecContext*)malloc(sizeof(AVCodecContext));
		//st_info->pCodec = (AVCodec*)malloc(sizeof(AVCodec));

		av_log_set_level(AV_LOG_QUIET);
		// Open video file
 		if(avformat_open_input(&st_info->pFormatCtx,st_info->filename,NULL,NULL)!=0){
			return -1 ; // Couldn't open file
		}
	 
		// Retrieve stream information
		if(avformat_find_stream_info(st_info->pFormatCtx,NULL)<0){
			return -1; // Couldn't find stream information
		}

		unsigned int i;

		// Find the video stream
		for(i=0; i< st_info->pFormatCtx->nb_streams; i++)
		{
			if(st_info->pFormatCtx->streams[i]->codec->codec_type==AVMEDIA_TYPE_VIDEO) 
			{
				st_info->videoStream=i;
				break;
			}
		}

		if(st_info->videoStream==-1){
		    return -1; //no video stream
		}
	
		// Get a pointer to the codec context for the video stream
		st_info->pCodecCtx = st_info->pFormatCtx->streams[st_info->videoStream]->codec;

		// Find the decoder
		st_info->pCodec = avcodec_find_decoder(st_info->pCodecCtx->codec_id);
		if(st_info->pCodec==NULL) 
		{
			return -1 ; // Codec not found
		}
		// Open codec
		if(avcodec_open2(st_info->pCodecCtx, st_info->pCodec,NULL)<0){
		    return -1; // Could not open codec
		}
		
		st_info->width = (st_info->width<=0) ? st_info->pCodecCtx->width : st_info->width;
		st_info->height = (st_info->height<=0) ? st_info->pCodecCtx->height : st_info->height;
		
	} 

	AVFrame *pFrame;

	// Allocate video frame
	pFrame=avcodec_alloc_frame();
		
	// Allocate an AVFrame structure
	AVFrame *pConvertedFrame = avcodec_alloc_frame();
	if(pConvertedFrame==NULL){
	  return -1;
	}
		
	uint8_t *buffer;
	int numBytes;
	// Determine required buffer size and allocate buffer
	numBytes=avpicture_get_size(ffmpeg_pixfmt, st_info->width,st_info->height);
	buffer=(uint8_t *)av_malloc(numBytes*sizeof(uint8_t));
	if (buffer == NULL){
	    return -1;
	}

	avpicture_fill((AVPicture *)pConvertedFrame,buffer,ffmpeg_pixfmt,st_info->width,st_info->height);
		
	int frameFinished;
	int size = 0;
	AVPacket packet;
	int result = 1;
	CImg<uint8_t> next_image;
	SwsContext *c = sws_getContext(st_info->pCodecCtx->width, st_info->pCodecCtx->height, st_info->pCodecCtx->pix_fmt, st_info->width, st_info->height, ffmpeg_pixfmt , SWS_BICUBIC, NULL, NULL, NULL);
	while ((result >= 0) && (size < st_info->nb_retrieval))
	{

		result = av_read_frame(st_info->pFormatCtx, &packet); 
		if (result < 0) 
			break;
		if(packet.stream_index == st_info->videoStream) {
			
		int channels = ffmpeg_pixfmt == PIX_FMT_GRAY8 ? 1 : 3;
 		AVPacket avpkt;
                av_init_packet(&avpkt);
                avpkt.data = packet.data;
                avpkt.size = packet.size;
                //
                // HACK for CorePNG to decode as normal PNG by default
                // same method used by ffmpeg
                avpkt.flags = AV_PKT_FLAG_KEY;

                avcodec_decode_video2(st_info->pCodecCtx, pFrame, &frameFinished,&avpkt);

		   // avcodec_decode_video(st_info->pCodecCtx, pFrame, &frameFinished,
		   //                      packet.data,packet.size);
 
		    if(frameFinished) {
		    	if (st_info->current_index == st_info->next_index)
		    	{
			    st_info->next_index += st_info->step;
			   
			    sws_scale(c, pFrame->data, pFrame->linesize, 0, st_info->pCodecCtx->height, pConvertedFrame->data, pConvertedFrame->linesize);
				   	
				next_image.assign(*pConvertedFrame->data, channels, st_info->width,st_info->height,1,true);
				next_image.permute_axes("yzcx");
				pFrameList->push_back(next_image);
				size++;
				   	 
		    	}    
				st_info->current_index++;
		    }
    	  }
    	        av_free_packet(&packet);
	}
	
	av_free(buffer);
	buffer = NULL;
	av_free(pConvertedFrame);
	pConvertedFrame = NULL;
	av_free(pFrame);
	pFrame = NULL;
	sws_freeContext(c);
	c = NULL;
	if (result < 0)
	{
		avcodec_close(st_info->pCodecCtx);
		avformat_close_input(&st_info->pFormatCtx);
		st_info->pCodecCtx = NULL;
		st_info->pFormatCtx = NULL;
		st_info->pCodec = NULL;
		st_info->width = -1;
		st_info->height = -1;
	}
	return size; 
}

int GetNumberStreams(const char *file)
{
	 AVFormatContext *pFormatCtx;
	 av_log_set_level(AV_LOG_QUIET);
	 av_register_all();
	// Open video file
	if (avformat_open_input(&pFormatCtx, file, NULL, NULL))
	  return -1 ; // Couldn't open file
		 
	// Retrieve stream information
	if(avformat_find_stream_info(pFormatCtx, NULL)<0)
	  return -1; // Couldn't find stream information
	int result = pFormatCtx->nb_streams;
	avformat_close_input(&pFormatCtx);
	return result;
}

long GetNumberVideoFrames(const char *file)
{
    long nb_frames = 0L;
	AVFormatContext *pFormatCtx;
    av_log_set_level(AV_LOG_QUIET);
	av_register_all();
	// Open video file
	if (avformat_open_input(&pFormatCtx, file, NULL, NULL))
	  return -1 ; // Couldn't open file
			 
	// Retrieve stream information
	if(avformat_find_stream_info(pFormatCtx, NULL)<0)
	  return -1; // Couldn't find stream information
		
	// Find the first video stream
	int videoStream=-1;
	for(unsigned int i=0; i<pFormatCtx->nb_streams; i++)
	{
	     if(pFormatCtx->streams[i]->codec->codec_type==AVMEDIA_TYPE_VIDEO) 
	     {
		    videoStream=i;
		    break;
             }
	}
	if(videoStream==-1)
	   return -1; // Didn't find a video stream
	AVStream *str = pFormatCtx->streams[videoStream];
		
        nb_frames = str->nb_frames;
	if (nb_frames > 0)
	{   //the easy way if value is already contained in struct 
	    avformat_close_input(&pFormatCtx);
	    return nb_frames;
	}
	else { // frames must be counted
	    AVPacket packet;
		nb_frames = (long)av_index_search_timestamp(str,str->duration, AVSEEK_FLAG_ANY|AVSEEK_FLAG_BACKWARD);
		// Close the video file
		 int timebase = str->time_base.den / str->time_base.num;
               if (nb_frames <= 0)
                       nb_frames = str->duration/timebase;
		avformat_close_input(&pFormatCtx); 
		return nb_frames;
	}
}

float fps(const char *filename)
{
        float result = 0;
	AVFormatContext *pFormatCtx;
	
	// Open video file
	if (avformat_open_input(&pFormatCtx, filename, NULL, NULL))
	  return -1 ; // Couldn't open file
				 
	// Retrieve stream information
	if(avformat_find_stream_info(pFormatCtx,NULL)<0)
	  return -1; // Couldn't find stream information
			
	// Find the first video stream
	int videoStream=-1;
	for(unsigned int i=0; i<pFormatCtx->nb_streams; i++)
	{
		     if(pFormatCtx->streams[i]->codec->codec_type==AVMEDIA_TYPE_VIDEO) 
		     {
			    videoStream=i;
			    break;
		     }
	}
        if(videoStream==-1)
	    return -1; // Didn't find a video stream
	
	int num = (pFormatCtx->streams[videoStream]->r_frame_rate).num;
	int den = (pFormatCtx->streams[videoStream]->r_frame_rate).den;
	result = num/den;

	avformat_close_input(&pFormatCtx);
	
	return result;

}
