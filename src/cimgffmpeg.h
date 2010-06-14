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



#ifndef CIMGFFMPEG_H_
#define CIMGFFMPEG_H_

#define cimg_display 0
#define cimg_debug 0

#include "CImg.h"

#define __STDC_CONSTANT_MACROS

extern "C" {
	#include "./libavformat/avformat.h"
	#include "./libavcodec/avcodec.h"
	#include "./libswscale/swscale.h"
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
