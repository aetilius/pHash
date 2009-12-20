#ifndef CIMGFFMPEG_H_
#define CIMGFFMPEG_H_

#include "CImg.h"

/*
   struct to hold video information to be passed to PlayVideo thread
*/
struct video_info_t {
	int fps;
	void *pList;
};

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

void vfinfo_close(VFInfo  *vfinfo){
    if (vfinfo->pFormatCtx != NULL){
	avcodec_close(vfinfo->pCodecCtx);
	vfinfo->pCodecCtx = NULL;
	av_close_input_file(vfinfo->pFormatCtx);
	vfinfo->pFormatCtx = NULL;
	vfinfo->width = -1;
	vfinfo->height = -1;
    }
}

int ReadFrames(VFInfo *st_info, CImgList<uint8_t> *pFrameList, unsigned int low_index, unsigned int hi_index);


int NextFrames(VFInfo *st_info, CImgList<uint8_t> *pFrameList);


int GetNumberStreams(const char *file);


long GetNumberVideoFrames(const char *file);

float fps(const char *filename);

#endif /*CIMGFFMPEG_H_*/
