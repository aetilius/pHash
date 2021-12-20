#ifndef CIMGFFMPEG_H_
#define CIMGFFMPEG_H_

#define cimg_display 0

#include "CImg.h"
#include "string.h"
#include "windows.h"

extern "C" {
#include "libavcodec/avcodec.h"
#include "libavformat/avformat.h"
#include "libswscale/swscale.h"
}
/*
   struct to hold video information to be passed to PlayVideo thread
*/
struct video_info_t {
    int fps;
    void *pList;
};

using namespace cimg_library;

/**
 *   !struct to hold information for reading from a media file
 *
 **/

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
    char *filename;
} VFInfo;

/** !close file when done getting frames from file **/
__declspec(dllexport) void vfinfo_close(VFInfo *vfinfo);

//! Read frames from the given file into the CImgList parameter
/*
 *  Reads nb_retrieval frames from the video stream from low_index to hi_index
 * in step number of intervals. The value of pixelformat parameter determines
 * the frame type. The pixel format value is translated into ffmpeg library's
 * corresponding PIX_FMT enumeration.
 *
 *  NOTE: opens and closes the file within function. No need to do so outside
 * function call.
 *
 *  @param st_info - state information for current file
 *  @param pFrameList - CImgList to append new frames into
 *  @param lo_index   int - the start index frame
 *  @param hi_index   int - the end index frame
 *  @param step       int - the interval
 *  @return int  -  the number of frames read; -1 for unable to read frames; 0
 * for no frames.
 */
__declspec(dllexport) int ReadFrames(VFInfo *st_info,
                                     CImgList<uint8_t> *pFrameList,
                                     unsigned int low_index,
                                     unsigned int hi_index);

//! read frames from stream at index where it left off on previous invocation of
//! function
/*
        Read up to nb_retrieval frames from filename video stream at step
   intervals; convert frames to given pixelformat; insert in empty pFrameList
   CImgList; meant to be called consecutively until the end of stream is
   reached, at which point another call will reinitialize at the beginning of
   the stream.  NOTE: Will open file at beginning and close file at end, no need
   to call vf_info_close() unless finished reading frames  before end of call.
        @param st_info - ptr to struct for state information.
        @param pFrameList CImgList<>* - an empty CImgList<> in which to append
   newly read frames
        @return int number of frames read or 0 if end of stream reached,  -1
   unable to read from stream
*/
__declspec(dllexport) int NextFrames(VFInfo *st_info,
                                     CImgList<uint8_t> *pFrameList);

//! Get number of streams contained in given file
/*
        Return the number of streams contained in the given file format.
        @param file const char*
        @return number of streams, -1 for error
*/
__declspec(dllexport) int GetNumberStreams(const char *file);

//! get number of video frames contained in file
/*
 * return number of video frames in file.
 * @param file const char* - video file
 * @return number of video frames, -1 for error
 */
__declspec(dllexport) long GetNumberVideoFrames(const char *file);

/** !determine frames per second
 *   @param  filename - name of file
 *   @return float value for frames per second, -1 for error
 **/
__declspec(dllexport) float fps(const char *filename);

/** !retrieve significant key frames from video
 *    @param filename - name of file
 *    @return CImgList - list of keyframes, NULL for error
 *
 **/
__declspec(dllexport) CImgList<uint8_t> *ph_getKeyFramesFromVideo(
    const char *filename);

#endif /*CIMGFFMPEG_H_*/
