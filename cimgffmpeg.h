/*
#	File: cimgffmpeg.h
#	
#	Description: plugin allowing convenient access to
#                 the frames of a video file through the
#                 CImg and CImgList classes; uses the ffmpeg
#                 function library.
#	Copyright: David G. Starkweather
#              starkdg@users.sourceforge.net
#              starkdg@comcast.net
#	License:
#   
#   Instructions:  This header depends on the ffmpeg libraries v13263
#
#		1.) obtain ffmpeg from svn (svn://svn.mplayerhq.hu/ffmpeg/trunk)
#
# 		2.) configure with the following options
#            ./configure --prefix=/usr/local --enable-shared --enable-swscale 
#                         --enable-avfilter --enable-pthreads --enable-gpl 
#                         --disable-ffmpeg --disable-ffserver --disable-ffplay
#           (I disable the ffmpeg, ffserver and ffplay utilities only because they
             are not needed by this header.)
#		3.) make
#  		4.) make install
#
#       Simply include this header in your code and link to the necessary libs.
#       You will probably have to edit the #include "cimgffmpeg.h" directives to
#       point to the correct directory level.
#
#       For your project that uses this header:
#
#       A command line invocation of the compiler might look something like this:
#
#		g++ -Dcimg_use_xshm -Dcimg_use_xrandr -Ucimg_use_xrandr -Ucimg_use_xshm -I/usr/include/X11 
#            -O2 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"src/ffmpegplugindemo.d" -MT"src/ffmpegplugindemo.d" -o"src/ffmpegplugindemo.o" "../src/ffmpegplugindemo.cpp"
#
#		And the linker:
#
#       g++ -L/usr/local/lib -L/usr/X11R6/lib64 -o{nameofyourproject}" ./src/{yourproject}.o
#             -lavutil -lavformat -lavcodec -lavutil -lswscale -lm -ldl -X11 -lXext -lXrandr 
#             -lpthread
#
#		set LD_LIBRARY_PATH environment variable to point to the lib location (e.g. /usr/local/lib)
*/

#ifndef CIMGFFMPEG_H_
#define CIMGFFMPEG_H_

#include <ctime>
#include "string.h"
#include <pthread.h>
#include "CImg.h"
extern "C" {
	#include "avformat.h"
	#include "avcodec.h"
	#include "swscale.h"
}
/*
   struct to hold video information to be passed to PlayVideo thread
*/
struct video_info_t {
	int fps;
	void *pList;
};

using namespace cimg_library;

//! Read frames from the given file into the CImgList parameter
/*  
 *  Reads nb_retrieval frames from the video stream from low_index to hi_index in step number of intervals.
 *  The value of pixelformat parameter determines the frame type. The pixel format value is translated into
 *  ffmpeg library's corresponding PIX_FMT enumeration.
 * 
 *  @param filename const char* - complete filename
 *  @param pFrameList CImgList<>* - an empty CImgList
 *  @param lo_index   int - the start index frame
 *  @param hi_index   int - the end index frame
 *  @param step       int - the interval 
 *  @param nb_retrieval int - maximum number of frames to retrieve
 *  @param pixelformat int - pixel format ( =0 for GRAY8 (8bpp grayscale) and =1 for RGB24 (24 bpp RGB)  )
 *  @return int  -  the number of frames read; -1 for unable to read frames; 0 for no frames.
*/
int ReadFrames(const char *filename, CImg<float> &videoframes, unsigned int low_index, unsigned int hi_index, int step = 1, long nb_retrieval=100, int pixelformat = 0)
{

	//determine destination pixel format from pixelformat input variable
	//ffmpeg_pixfmt values are taken from enum PixelFormat in avutil.h
	int ffmpeg_pixfmt;
	if (pixelformat == 0)
		ffmpeg_pixfmt = PIX_FMT_RGB24;
	else //(pixelformat == 1)
		ffmpeg_pixfmt = PIX_FMT_GRAY8;
	
	av_register_all();
	AVFormatContext *pFormatCtx;

        int index_vf = 0;	

	// Open video file
	if(av_open_input_file(&pFormatCtx, filename, NULL, 0, NULL)!=0)
	  return -1 ; // Couldn't open file
	 
	// Retrieve stream information
	if(av_find_stream_info(pFormatCtx)<0)
	  return -1; // Couldn't find stream information
	
	//dump_format(pFormatCtx,0,NULL,0);//debugging function to print infomation about format
	
	unsigned int i;
	AVCodecContext *pCodecCtx;

	// Find the video stream
	int videoStream=-1;
	for(i=0; i<pFormatCtx->nb_streams; i++)
	{
	     if(pFormatCtx->streams[i]->codec->codec_type==CODEC_TYPE_VIDEO) 
	     {
		    videoStream=i;
		    break;
		  }
		if(videoStream==-1)
		  return -1; //no video stream
	}
	
	// Get a pointer to the codec context for the video stream
	pCodecCtx=pFormatCtx->streams[videoStream]->codec;

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
		
	AVFrame *pFrame;

	// Allocate video frame
	pFrame=avcodec_alloc_frame();
		
	// Allocate an AVFrame structure
	AVFrame *pConvertedFrame = avcodec_alloc_frame();
	if(pConvertedFrame==NULL)
	  return -1;
		
	uint8_t *buffer;
	int numBytes;
	// Determine required buffer size and allocate buffer
	numBytes=avpicture_get_size(ffmpeg_pixfmt, pCodecCtx->width,pCodecCtx->height);
	buffer=(uint8_t *)av_malloc(numBytes*sizeof(uint8_t));

	avpicture_fill((AVPicture *)pConvertedFrame, buffer, ffmpeg_pixfmt,pCodecCtx->width, pCodecCtx->height);
		
	int frameFinished;
	unsigned int size = 0;
	unsigned int current_index = low_index;
	unsigned int next_index = current_index;
	AVPacket packet;
	
	while(av_read_frame(pFormatCtx, &packet)>=0) 
	{
    	  if(packet.stream_index==videoStream) {
			
		    avcodec_decode_video(pCodecCtx, pFrame, &frameFinished,
		                         packet.data, packet.size);
		    
		    if(frameFinished) {
		    	 
		    	if (current_index == next_index)
		    	{
		    		
		    		next_index += step;
    	    		SwsContext *c = sws_getContext(pCodecCtx->width, pCodecCtx->height, pCodecCtx->pix_fmt, pCodecCtx->width, pCodecCtx->height, ffmpeg_pixfmt , 1, NULL, NULL, NULL);
				   	sws_scale(c, pFrame->data, pFrame->linesize, 0, pCodecCtx->height, pConvertedFrame->data, pConvertedFrame->linesize);
				   	
				   	if (ffmpeg_pixfmt == PIX_FMT_RGB24)
				   	{
				   		CImg<unsigned char> *pNextImage = new CImg<unsigned char>(*pConvertedFrame->data,3,pCodecCtx->width,pCodecCtx->height,1,true);
				   		CImg<unsigned char> NextImage = pNextImage->get_permute_axes("yzvx");
				   	}
				   	else if (ffmpeg_pixfmt == PIX_FMT_GRAY8){
				   		CImg<uint8_t> *pNextImage = new CImg<uint8_t>(*pConvertedFrame->data,1,pCodecCtx->width,pCodecCtx->height,1,true);
				   		CImg<float> NextImage = pNextImage->get_permute_axes("yzvx");
						NextImage.blur(1.0).resize(32,32);
				   		
						cimg_forXY(NextImage,X,Y){
						    videoframes(X,Y,index_vf) = NextImage(X,Y);
						}
						index_vf++;
						
				   	}
				   	size++;
		    	}    
				current_index++;
		    }
		    
		  av_free_packet(&packet);
		  if (next_index >= hi_index)
			  break;
    	  }
	}
	av_free(buffer);
	av_free(pConvertedFrame);
	av_free(pFrame);
	avcodec_close(pCodecCtx);
	av_close_input_file(pFormatCtx);
	return size; 
}
//! read frames from stream at index where it left off on previous invocation of function
/*
 	Read up to nb_retrieval frames from filename video stream at step intervals; convert 
 	frames to given pixelformat; insert in empty pFrameList CImgList; meant to be called consecutively
 	until the end of stream is reached, at which point another call will reinitialize at the beginning
 	of the stream.
 	@param filename const char* - name of the file
 	@param pFrameList CImgList<>* - an empty CImgList<> 
 	@param step int - the step interval
 	@param nb_retrieval int - max number to retrieve at a time
 	@param pixelformat int - pixel format ( 0 = GRAY8 (8bpp grayscale), 1 = RGB24 (24 bpp RGB) )
 	@return number of frames read or 0 if end of stream reached,  -1 unable to read from stream
*/
int NextFrames(const char *filename, CImgList<unsigned char> *pFrameList, int step = 1, unsigned int nb_retrieval=100, int pixelformat = 0)
{
	static int initd = 0;
	static int videoStream = -1;
	static AVFormatContext *pFormatCtx = (AVFormatContext*)malloc(sizeof(AVFormatContext));
	static AVCodecContext *pCodecCtx = (AVCodecContext*)malloc(sizeof(AVCodecContext));
	static AVCodec *pCodec = (AVCodec*)malloc(sizeof(AVCodec));
	
	//determine destination pixel format from pixelformat input variable
	//ffmpeg_pixfmt values are taken from enum PixelFormat in avutil.h
	int ffmpeg_pixfmt = PIX_FMT_RGB24;
	if (pixelformat == 1)
		ffmpeg_pixfmt = PIX_FMT_GRAY8;
	
	if (!initd)
	{
		av_register_all();
		
		// Open video file
		if(av_open_input_file(&pFormatCtx, filename, NULL, 0, NULL)!=0)
			return -1 ; // Couldn't open file
	 
		// Retrieve stream information
		if(av_find_stream_info(pFormatCtx)<0)
			return -1; // Couldn't find stream information
	
		unsigned int i;

		// Find the video stream
		for(i=0; i<pFormatCtx->nb_streams; i++)
		{
			if(pFormatCtx->streams[i]->codec->codec_type==CODEC_TYPE_VIDEO) 
			{
				videoStream=i;
				break;
			}
		if(videoStream==-1)
		  return -1; //no video stream
		}
	
		// Get a pointer to the codec context for the video stream
		pCodecCtx=pFormatCtx->streams[videoStream]->codec;

		// Find the decoder
		pCodec=avcodec_find_decoder(pCodecCtx->codec_id);
		if(pCodec==NULL) 
		{
			return -1 ; // Codec not found
		}
		// Open codec
		if(avcodec_open(pCodecCtx, pCodec)<0)
			return -1; // Could not open codec
		
		initd = 1;
	}
	AVFrame *pFrame;

	// Allocate video frame
	pFrame=avcodec_alloc_frame();
		
	// Allocate an AVFrame structure
	AVFrame *pConvertedFrame = avcodec_alloc_frame();
	if(pConvertedFrame==NULL)
	  return -1;
		
	uint8_t *buffer;
	int numBytes;
	// Determine required buffer size and allocate buffer
	numBytes=avpicture_get_size(ffmpeg_pixfmt, pCodecCtx->width,pCodecCtx->height);
	buffer=(uint8_t *)av_malloc(numBytes*sizeof(uint8_t));

	avpicture_fill((AVPicture *)pConvertedFrame, buffer, ffmpeg_pixfmt,pCodecCtx->width, pCodecCtx->height);
		
	int frameFinished;
	unsigned int size = 0;
	unsigned int current_index = 0;
	unsigned int next_index = current_index;
	AVPacket packet;
	int result = 1;
	while ((result >= 0) && (size < nb_retrieval))
	{
		result = av_read_frame(pFormatCtx, &packet); 
		if (result < 0) 
			break;
		if(packet.stream_index==videoStream) {
			
		    avcodec_decode_video(pCodecCtx, pFrame, &frameFinished,
		                         packet.data, packet.size);
 
		    if(frameFinished) {
		    	if (current_index == next_index)
		    	{
		    		
		    		next_index += step;
    	    		SwsContext *c = sws_getContext(pCodecCtx->width, pCodecCtx->height, pCodecCtx->pix_fmt, pCodecCtx->width, pCodecCtx->height, ffmpeg_pixfmt , 1, NULL, NULL, NULL);
				   	sws_scale(c, pFrame->data, pFrame->linesize, 0, pCodecCtx->height, pConvertedFrame->data, pConvertedFrame->linesize);
				   	
				   	if (ffmpeg_pixfmt == PIX_FMT_RGB24)
				   	{   
				   		CImg<unsigned char> *pNextImage = new CImg<unsigned char>(*pConvertedFrame->data,3,pCodecCtx->width,pCodecCtx->height,1,true);
				   		CImg<unsigned char> NextImage = pNextImage->get_permute_axes("yzvx");
				   		pFrameList->push_back(NextImage);
				   		size++;
				   	}
				   	else if (ffmpeg_pixfmt == PIX_FMT_GRAY8){
				   		
				   		CImg<unsigned char> *pNextImage = new CImg<unsigned char>(*pConvertedFrame->data,1,pCodecCtx->width,pCodecCtx->height,1,true);
				   		CImg<unsigned char> NextImage = pNextImage->get_permute_axes("yzvx");
				   		pFrameList->push_back(NextImage);
				   		size++;
				   	}
				   	 
		    	}    
				current_index++;
		    }
    	  }
    	  
		  av_free_packet(&packet);
	}
	
	 
	av_free(buffer);
	av_free(pConvertedFrame);
	av_free(pFrame);
	if (result < 0)
	{
		avcodec_close(pCodecCtx);
		av_close_input_file(pFormatCtx);
		initd = 0;
		return 0;
	}
	return size; 
}

//! Get number of streams contained in given file
/*
	Return the number of streams contained in the given file format.
	@param file const char*
	@return number of streams
*/
int GetNumberStreams(const char *file)
{
	 AVFormatContext *pFormatCtx;
	 av_register_all();
	// Open video file
	if (av_open_input_file(&pFormatCtx, file, NULL, 0, NULL))
	  return -1 ; // Couldn't open file
		 
	// Retrieve stream information
	if(av_find_stream_info(pFormatCtx)<0)
	  return -1; // Couldn't find stream information
	int result = pFormatCtx->nb_streams;
	av_close_input_file(pFormatCtx);
	return result;
}
//! get number of video frames contained in file
/*
 * return number of video frames in file.
 * @param file const char* - video file
 * @return number of video frames
*/
long GetNumberVideoFrames(const char *file)
{
        long nb_frames = 0L;
	AVFormatContext *pFormatCtx;
	av_register_all();
	// Open video file
	if (av_open_input_file(&pFormatCtx, file, NULL, 0, NULL))
	  return -1 ; // Couldn't open file
			 
	// Retrieve stream information
	if(av_find_stream_info(pFormatCtx)<0)
	  return -1; // Couldn't find stream information
		
	// Find the first video stream
	int videoStream=-1;
	for(unsigned int i=0; i<pFormatCtx->nb_streams; i++)
	{
	     if(pFormatCtx->streams[i]->codec->codec_type==CODEC_TYPE_VIDEO) 
	     {
		    videoStream=i;
		    break;
             }
	     if(videoStream==-1)
		   return -1; // Didn't find a video stream
	}
		
        nb_frames = pFormatCtx->streams[videoStream]->nb_frames;
	if (nb_frames > 0)
	{   //the easy way if value is already contained in struct 
	    av_close_input_file(pFormatCtx);
	    return nb_frames;
	}
	else { // frames must be counted
	        AVPacket packet;
		
		//read each frame - one frame per packet
		while(av_read_frame(pFormatCtx, &packet)>=0) 
		{
		  if(packet.stream_index==videoStream) {
                  //packet is from video stream
		      nb_frames++;
		  }
		  av_free_packet(&packet);
		}

		// Close the video file
		av_close_input_file(pFormatCtx); 

		return nb_frames;
	}
}
float fps(const char *filename)
{
	AVFormatContext *pFormatCtx;
	
	// Open video file
	if (av_open_input_file(&pFormatCtx, filename, NULL, 0, NULL))
	  return -1 ; // Couldn't open file
				 
	// Retrieve stream information
	if(av_find_stream_info(pFormatCtx)<0)
	  return -1; // Couldn't find stream information
			
	// Find the first video stream
	int videoStream=-1;
	for(unsigned int i=0; i<pFormatCtx->nb_streams; i++)
	{
		     if(pFormatCtx->streams[i]->codec->codec_type==CODEC_TYPE_VIDEO) 
		     {
			    videoStream=i;
			    break;
			  }
		     if(videoStream==-1)
			   return -1; // Didn't find a video stream
	}
	av_close_input_file(pFormatCtx);		
	
	int num = (pFormatCtx->streams[videoStream]->r_frame_rate).num;
	int den = (pFormatCtx->streams[videoStream]->r_frame_rate).den;
	
	
	return (num/den);

}

//! play a list of CImg's given a CImgList.   
  
/*
 *  This function is not really meant to function as a video player as much as a simple 
 *  debugging tool.  It is meant to be run in a thread. The video is played at 1 CImg per second
 *  @param pVideo - a pointer to a CImgList object containing the CImg's to be played 
*/
void *PlayVideo(void *info)
{
	video_info_t *pvideo_info = (video_info_t *)info;
	int fps = pvideo_info->fps;
	CImgList<unsigned char> *videolist = (CImgList<unsigned char> *)pvideo_info->pList;
	CImg<> current = videolist->front();
    
	CImgDisplay disp(current,"video");
	
	//display frame every second
	clock_t interval = (time_t)((1/fps)*(CLOCKS_PER_SEC)); 
			 
	clock_t n = clock();
	while (!videolist->is_empty() && !disp.is_closed)
	{
		if (videolist->size > 1)
		{
			videolist->pop_front();
			current = videolist->front();
			disp.display(current);
			while ( n + interval > clock()){};
			n = clock();
		}
		else if (videolist->size == 1){
			videolist->pop_front();
		}
		else
			pthread_exit(NULL);
	}
	pthread_exit(NULL);
}
/**
 * !write a sequence of CImg's to video file.
 * Purpose: write sequence of images  in cimg list at a given framerate 
 * using given filename, use pixel format PIX_FMT_YUV420P
 * The function expects the images to be either gray scale (dim v = 1) or 
 * rgb images (dim v == 3).  Other sizes return -1.
 * 
 * @param pImlist reference to CImgList
 * @param filename name of video file to create - with proper extension - e.g. 
 *                                ".avi, .mpeg, .wmv, etc."
 * @param fps int value for number of frames per second
 * @param stream_duration maximum length of video stream (secs)
 * @return int number of frames recorded to file, < 0 for errors
 **/
template<class T>
int WriteFrames(CImgList<T> *pImList,const char *filename, int fps,float stream_duration=100.0)
{
    int width  = (pImList->front()).dimx();
    int height = (pImList->front()).dimy();
    int nb_channels = (pImList->front()).dimv();
    if ((nb_channels != 1) && (nb_channels != 3))
	throw new CImgIOException("dimv of image not acceptable");
    if (!filename)
        throw new CImgIOException("no file name given");
    int frame_count = 0;
    PixelFormat dest_pxl_fmt = PIX_FMT_YUV420P;
    PixelFormat src_pxl_fmt  = (nb_channels == 3) ? PIX_FMT_RGB24 : PIX_FMT_GRAY8;
    
    int sws_flags = SWS_FAST_BILINEAR;//interpolation method (keeping same size images for now)

    AVOutputFormat *fmt = NULL;

    fmt = guess_format(NULL,filename,NULL);
    if (!fmt){
	//default format "mpeg" 
        fmt = guess_format("mpeg",NULL,NULL);
    }
    if (!fmt){
	//unable to retrieve output format
        throw new CImgIOException("could not determine format from filename");
    }

    AVFormatContext *oc = NULL;
    oc = av_alloc_format_context(); 
    if (!oc){
	//unable to allocate format context
        throw new CImgIOException("mem allocation error for format context");
    }

    AVCodec *codec = NULL;
    AVFrame *picture = NULL;
    AVFrame *tmp_pict = NULL;
    oc->oformat = fmt;
    snprintf(oc->filename, sizeof(oc->filename),"%s",filename);

    av_register_all();

    //add video stream
    int stream_index = 0;
    AVStream *video_str = NULL;
    if (fmt->video_codec != CODEC_ID_NONE) {
	video_str = av_new_stream(oc,stream_index);
	if (!video_str){
	    //no stream allocated
            av_free(oc);
	    throw new CImgIOException("unable to create new video stream");
	}
    } else {
        //no codec identified
        av_free(oc);
	throw new CImgIOException("no codec identified"); 
    }

    AVCodecContext *c = video_str->codec;
    c->codec_id = fmt->video_codec;
    c->codec_type = CODEC_TYPE_VIDEO;
    c->bit_rate = 400000;
    c->width = width;
    c->height = height;
    c->time_base.num = 1;
    c->time_base.den = fps;
    c->gop_size = 12;
    c->pix_fmt = dest_pxl_fmt;
    if (c->codec_id == CODEC_ID_MPEG2VIDEO)
	c->max_b_frames = 2;
    if (c->codec_id == CODEC_ID_MPEG1VIDEO)
        c->mb_decision = 2;
 
    if (av_set_parameters(oc,NULL) < 0){
        //parameters not properly set
        av_free(oc);
        throw new CImgIOException("parameters for avcodec not properly set");
    }
   
    //dump_format(oc,0,filename,1);

    //open codecs and alloc buffers
    codec = avcodec_find_encoder(c->codec_id);
    if (!codec){
	//unable to find codec
        av_free(oc);
        throw new CImgIOException("no codec found");
    }
    if (avcodec_open(c,codec) < 0)
        //fail to open codec
	throw new CImgIOException("unable to open codec");

    tmp_pict = avcodec_alloc_frame();
    if (!tmp_pict){
	//unable to allocate memory for tmp_pict frame
        avcodec_close(video_str->codec);
        av_free(oc);
        throw new CImgIOException("mem alloc error for tmp_pict data buffer");
    }
    tmp_pict->linesize[0] = (src_pxl_fmt == PIX_FMT_RGB24) ? 3*width : width ;
    tmp_pict->type = FF_BUFFER_TYPE_USER;
    int	tmp_size = avpicture_get_size(src_pxl_fmt,width,height);
    uint8_t *tmp_buffer = (uint8_t*)av_malloc(tmp_size);
    if (!tmp_buffer){
	//unable to allocate memory for tmp buffer
        av_free(tmp_pict);
        avcodec_close(video_str->codec);
        av_free(oc);
        throw new CImgIOException("mem alloc error for tmp buffer");
    }

    //associate buffer with tmp_pict
    avpicture_fill((AVPicture*)tmp_pict,tmp_buffer,src_pxl_fmt,width,height);

    picture = avcodec_alloc_frame();
    if (!picture){
	//unable to allocate picture frame
        av_free(tmp_pict->data[0]);
        av_free(tmp_pict);
        avcodec_close(video_str->codec);
        av_free(oc);
        throw new CImgIOException("mem alloc error for picture frame");
    }

    int size = avpicture_get_size(c->pix_fmt,width,height);
    uint8_t *buffer = (uint8_t*)av_malloc(size);
    if (!buffer){
	//unable to allocate buffer
        av_free(picture);
        av_free(tmp_pict->data[0]);
        av_free(tmp_pict);
        avcodec_close(video_str->codec);
        av_free(oc);
        throw new CImgIOException("mem alloc error for picture frame buffer");
    }
    //associate the buffer with picture
    avpicture_fill((AVPicture*)picture,buffer,c->pix_fmt,width,height);

    //open file
    if ( !(fmt->flags & AVFMT_NOFILE) ){
	if (url_fopen(&oc->pb,filename,URL_WRONLY) < 0)
	    throw new CImgIOException("unable to open file");
    }

    if (av_write_header(oc) < 0)
	throw new CImgIOException("could not write header");
    
    double video_pts;

    SwsContext *img_convert_context = NULL;
    img_convert_context = sws_getContext(width,height,src_pxl_fmt,
                                         c->width,c->height,c->pix_fmt,sws_flags,NULL,NULL,NULL);
    if (!img_convert_context){
	//unable to get swscale context
        if (!(fmt->flags & AVFMT_NOFILE))
	    url_fclose(oc->pb);
        av_free(picture->data);
        av_free(picture);
        av_free(tmp_pict->data[0]);
        av_free(tmp_pict);
        avcodec_close(video_str->codec);
        av_free(oc);
        throw new CImgIOException("unable to get conversion context");
    }
    int ret=0, out_size;
    uint8_t *video_outbuf = NULL;
    int video_outbuf_size = 1000000;
    video_outbuf = (uint8_t*)av_malloc(video_outbuf_size);    
    if (!video_outbuf){
	if (!(fmt->flags & AVFMT_NOFILE))
	    url_fclose(oc->pb);
	av_free(picture->data);
        av_free(picture);
        av_free(tmp_pict->data[0]);
        av_free(tmp_pict);
        avcodec_close(video_str->codec);
        av_free(oc);
        throw new CImgIOException("mem alloc error");
    }

    //loop through each image in list
    for (unsigned int i = 0; i < pImList->size; i++)
    {
	frame_count++;
        CImg<uint8_t> currentIm = pImList->at(i);
	CImg<uint8_t> red,green,blue,gray;
	if (src_pxl_fmt == PIX_FMT_RGB24){
	    red = currentIm.get_channel(0);
	    green = currentIm.get_channel(1);
	    blue = currentIm.get_channel(2);
            cimg_forXY(red,X,Y){
		//assign pizel values to data buffer in interlaced RGBRGB ... format
		tmp_pict->data[0][Y*tmp_pict->linesize[0] + 3*X]     = red(X,Y);
                tmp_pict->data[0][Y*tmp_pict->linesize[0] + 3*X + 1] = green(X,Y);
                tmp_pict->data[0][Y*tmp_pict->linesize[0] + 3*X + 2] = blue(X,Y);
	    }
	} else {
	    gray = currentIm.get_channel(0);
	    cimg_forXY(gray,X,Y){
                tmp_pict->data[0][Y*tmp_pict->linesize[0] + X] = gray(X,Y);
	    }
	}


	if (video_str)
	    video_pts = (video_str->pts.val * video_str->time_base.num)/(video_str->time_base.den);
        else
	    video_pts = 0.0;
        if ((!video_str) || (video_pts >= stream_duration))
	    break;

	if (sws_scale(img_convert_context,tmp_pict->data,tmp_pict->linesize,0,c->height,picture->data,picture->linesize) < 0){
            break;;//break out of loop
	}
	out_size = avcodec_encode_video(c,video_outbuf,video_outbuf_size,picture);
        
        if (out_size > 0){
	    AVPacket pkt;
            av_init_packet(&pkt);
            pkt.pts = av_rescale_q(c->coded_frame->pts,c->time_base,video_str->time_base);
            if (c->coded_frame->key_frame){
		pkt.flags |= PKT_FLAG_KEY;
	    }
            pkt.stream_index = video_str->index;
            pkt.data = video_outbuf;
            pkt.size = out_size;
            ret = av_write_frame(oc,&pkt);
	} else if (out_size < 0){
            break;//failure occured in avcodec_encode_video() 
	}
        if (ret != 0){
            //error occured in writing frame
            break;
	}
 
    }
    //close codec
    if (video_str){
	avcodec_close(video_str->codec);
        av_free(picture->data[0]);
        av_free(picture);
        av_free(tmp_pict->data[0]);
        av_free(tmp_pict);
    }
    if (av_write_trailer(oc) < 0)
	throw new CImgIOException("could not write trailer");
    av_freep(&oc->streams[stream_index]->codec);
    av_freep(&oc->streams[stream_index]);
    if (!(fmt->flags & AVFMT_NOFILE)){
        if (url_fclose(oc->pb) < 0)
	    throw new CImgIOException("could not close file");
    }
    av_free(oc);
    av_free(video_outbuf);

    return frame_count;
}
#endif /*CIMGFFMPEG_H_*/
