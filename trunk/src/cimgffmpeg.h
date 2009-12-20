#ifndef CIMGFFMPEG_H_
#define CIMGFFMPEG_H_

#define cimg_display 0
#define cimg_debug 0

#include "CImg.h"

extern "C" {
	#include "libavformat/avformat.h"
	#include "libavcodec/avcodec.h"
	#include "libswscale/swscale.h"
}

using namespace cimg_library;

typedef struct vf_info {
    int step;
    int nb_retrieval;
    int pixelformat;
    int videoStream;
    int width, height;
    long current_index;
    long next_index;
    AVFormatContext *pFormatCtx;
    AVCodecContext *pCodecCtx;
    AVCodec *pCodec;
    const char *filename;
} VFInfo;

void vfinfo_close(VFInfo  *vfinfo);

int ReadFrames(VFInfo *st_info, CImgList<uint8_t> *pFrameList, unsigned int low_index, unsigned int hi_index);


int NextFrames(VFInfo *st_info, CImgList<uint8_t> *pFrameList);


int GetNumberStreams(const char *file);


long GetNumberVideoFrames(const char *file);

float fps(const char *filename);

#endif /*CIMGFFMPEG_H_*/
