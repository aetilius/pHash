#include "cimgffmpeg.h"

#include <cmath>

__declspec(dllexport) void vfinfo_close(VFInfo *vfinfo) {
    if (vfinfo->pFormatCtx != NULL) {
        avcodec_close(vfinfo->pCodecCtx);
        vfinfo->pCodecCtx = NULL;
        av_close_input_file(vfinfo->pFormatCtx);
        vfinfo->pFormatCtx = NULL;
        vfinfo->width = -1;
        vfinfo->height = -1;
    }
}

_declspec(dllexport) int ReadFrames(VFInfo *st_info,
                                    CImgList<uint8_t> *pFrameList,
                                    unsigned int low_index,
                                    unsigned int hi_index) {
    // target pixel format
    PixelFormat ffmpeg_pixfmt;
    if (st_info->pixelformat == 0)
        ffmpeg_pixfmt = PIX_FMT_GRAY8;
    else
        ffmpeg_pixfmt = PIX_FMT_RGB24;

    st_info->next_index = low_index;

    if (st_info->pFormatCtx == NULL) {
        st_info->current_index = 0;

        av_register_all();
        av_log_set_level(AV_LOG_DEBUG);
        // Open video file
        if (av_open_input_file(&st_info->pFormatCtx, st_info->filename, NULL, 0,
                               NULL) != 0)
            return -1;  // Couldn't open file

        // Retrieve stream information
        if (av_find_stream_info(st_info->pFormatCtx) < 0)
            return -1;  // Couldn't find stream information

        // dump_format(pFormatCtx,0,NULL,0);//debugging function to print
        // infomation about format

        unsigned int i;
        // Find the video stream
        for (i = 0; i < st_info->pFormatCtx->nb_streams; i++) {
            if (st_info->pFormatCtx->streams[i]->codec->codec_type ==
                CODEC_TYPE_VIDEO) {
                st_info->videoStream = i;
                break;
            }
        }
        if (st_info->videoStream == -1) return -1;  // no video stream

        // Get a pointer to the codec context for the video stream
        st_info->pCodecCtx =
            st_info->pFormatCtx->streams[st_info->videoStream]->codec;
        if (st_info->pCodecCtx == NULL) {
            return -1;
        }

        // Find the decoder
        st_info->pCodec = avcodec_find_decoder(st_info->pCodecCtx->codec_id);
        if (st_info->pCodec == NULL) {
            return -1;  // Codec not found
        }
        // Open codec
        if (avcodec_open(st_info->pCodecCtx, st_info->pCodec) < 0)
            return -1;  // Could not open codec

        st_info->height = (st_info->height <= 0) ? st_info->pCodecCtx->height
                                                 : st_info->height;
        st_info->width =
            (st_info->width <= 0) ? st_info->pCodecCtx->width : st_info->width;
    }

    AVFrame *pFrame;

    // Allocate video frame
    pFrame = avcodec_alloc_frame();
    if (pFrame == NULL) return -1;

    // Allocate an AVFrame structure
    AVFrame *pConvertedFrame = avcodec_alloc_frame();
    if (pConvertedFrame == NULL) return -1;

    uint8_t *buffer;
    int numBytes;
    // Determine required buffer size and allocate buffer
    numBytes =
        avpicture_get_size(ffmpeg_pixfmt, st_info->width, st_info->height);
    buffer = (uint8_t *)av_malloc(numBytes * sizeof(uint8_t));
    if (buffer == NULL) return -1;

    avpicture_fill((AVPicture *)pConvertedFrame, buffer, ffmpeg_pixfmt,
                   st_info->width, st_info->height);

    int frameFinished;
    int size = 0;

    AVPacket packet;
    int result = 1;
    CImg<uint8_t> next_image;
    // CImg<uint8_t> tmp_img;
    SwsContext *c = sws_getContext(
        st_info->pCodecCtx->width, st_info->pCodecCtx->height,
        st_info->pCodecCtx->pix_fmt, st_info->width, st_info->height,
        ffmpeg_pixfmt, SWS_BICUBIC, NULL, NULL, NULL);
    while ((result >= 0) && (size < st_info->nb_retrieval) &&
           (st_info->current_index <= ((int)hi_index))) {
        result = av_read_frame(st_info->pFormatCtx, &packet);
        if (result < 0) break;
        if (packet.stream_index == st_info->videoStream) {
            avcodec_decode_video(st_info->pCodecCtx, pFrame, &frameFinished,
                                 packet.data, packet.size);
            if (frameFinished) {
                if (st_info->current_index == st_info->next_index) {
                    st_info->next_index += st_info->step;
                    sws_scale(c, pFrame->data, pFrame->linesize, 0,
                              st_info->pCodecCtx->height, pConvertedFrame->data,
                              pConvertedFrame->linesize);

                    if (ffmpeg_pixfmt == PIX_FMT_GRAY8) {
                        next_image.assign(*pConvertedFrame->data, 1,
                                          st_info->width, st_info->height, 1,
                                          true);
                        next_image.permute_axes("yzcx");
                        pFrameList->push_back(next_image);
                        size++;
                    } else if (ffmpeg_pixfmt == PIX_FMT_RGB24) {
                        next_image.assign(*pConvertedFrame->data, 3,
                                          st_info->width, st_info->height, 1,
                                          true);
                        next_image.permute_axes("yzcx");
                        pFrameList->push_back(next_image);
                        size++;
                    }
                }
                st_info->current_index++;
            }
            av_free_packet(&packet);
        }
    }

    if (result < 0) {
        avcodec_close(st_info->pCodecCtx);
        av_close_input_file(st_info->pFormatCtx);
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

__declspec(dllexport) int NextFrames(VFInfo *st_info,
                                     CImgList<uint8_t> *pFrameList) {
    PixelFormat ffmpeg_pixfmt;
    if (st_info->pixelformat == 0)
        ffmpeg_pixfmt = PIX_FMT_GRAY8;
    else
        ffmpeg_pixfmt = PIX_FMT_RGB24;

    if (st_info->pFormatCtx == NULL) {
        st_info->current_index = 0;
        st_info->next_index = 0;
        av_register_all();
        st_info->videoStream = -1;
        // st_info->pFormatCtx =
        // (AVFormatContext*)malloc(sizeof(AVFormatContext)); st_info->pCodecCtx
        // = (AVCodecContext*)malloc(sizeof(AVCodecContext)); st_info->pCodec =
        // (AVCodec*)malloc(sizeof(AVCodec));

        av_log_set_level(AV_LOG_QUIET);
        // Open video file
        if (av_open_input_file(&(st_info->pFormatCtx), st_info->filename, NULL,
                               0, NULL) != 0) {
            return -1;  // Couldn't open file
        }

        // Retrieve stream information
        if (av_find_stream_info(st_info->pFormatCtx) < 0) {
            return -1;  // Couldn't find stream information
        }

        unsigned int i;

        // Find the video stream
        for (i = 0; i < st_info->pFormatCtx->nb_streams; i++) {
            if (st_info->pFormatCtx->streams[i]->codec->codec_type ==
                CODEC_TYPE_VIDEO) {
                st_info->videoStream = i;
                break;
            }
        }

        if (st_info->videoStream == -1) {
            return -1;  // no video stream
        }

        // Get a pointer to the codec context for the video stream
        st_info->pCodecCtx =
            st_info->pFormatCtx->streams[st_info->videoStream]->codec;

        // Find the decoder
        st_info->pCodec = avcodec_find_decoder(st_info->pCodecCtx->codec_id);
        if (st_info->pCodec == NULL) {
            return -1;  // Codec not found
        }
        // Open codec
        if (avcodec_open(st_info->pCodecCtx, st_info->pCodec) < 0) {
            return -1;  // Could not open codec
        }

        st_info->width =
            (st_info->width <= 0) ? st_info->pCodecCtx->width : st_info->width;
        st_info->height = (st_info->height <= 0) ? st_info->pCodecCtx->height
                                                 : st_info->height;
    }

    AVFrame *pFrame;

    // Allocate video frame
    pFrame = avcodec_alloc_frame();

    // Allocate an AVFrame structure
    AVFrame *pConvertedFrame = avcodec_alloc_frame();
    if (pConvertedFrame == NULL) {
        return -1;
    }

    uint8_t *buffer;
    int numBytes;
    // Determine required buffer size and allocate buffer
    numBytes =
        avpicture_get_size(ffmpeg_pixfmt, st_info->width, st_info->height);
    buffer = (uint8_t *)av_malloc(numBytes * sizeof(uint8_t));
    if (buffer == NULL) {
        return -1;
    }

    avpicture_fill((AVPicture *)pConvertedFrame, buffer, ffmpeg_pixfmt,
                   st_info->width, st_info->height);

    int frameFinished;
    int size = 0;
    AVPacket packet;
    int result = 1;
    CImg<uint8_t> next_image;
    SwsContext *c = sws_getContext(
        st_info->pCodecCtx->width, st_info->pCodecCtx->height,
        st_info->pCodecCtx->pix_fmt, st_info->width, st_info->height,
        ffmpeg_pixfmt, SWS_BICUBIC, NULL, NULL, NULL);
    while ((result >= 0) && (size < st_info->nb_retrieval)) {
        result = av_read_frame(st_info->pFormatCtx, &packet);
        if (result < 0) break;
        if (packet.stream_index == st_info->videoStream) {
            avcodec_decode_video(st_info->pCodecCtx, pFrame, &frameFinished,
                                 packet.data, packet.size);

            if (frameFinished) {
                if (st_info->current_index == st_info->next_index) {
                    st_info->next_index += st_info->step;

                    sws_scale(c, pFrame->data, pFrame->linesize, 0,
                              st_info->pCodecCtx->height, pConvertedFrame->data,
                              pConvertedFrame->linesize);

                    if (ffmpeg_pixfmt == PIX_FMT_RGB24) {
                        next_image.assign(pConvertedFrame->data[0], 3,
                                          st_info->width, st_info->height, 1,
                                          true);
                        next_image.permute_axes("yzcx");
                        pFrameList->push_back(next_image);
                        size++;
                    } else if (ffmpeg_pixfmt == PIX_FMT_GRAY8) {
                        next_image.assign(pConvertedFrame->data[0], 1,
                                          st_info->width, st_info->height, 1,
                                          true);
                        next_image.permute_axes("yzcx");
                        pFrameList->push_back(next_image);
                        size++;
                    }
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
    if (result < 0) {
        avcodec_close(st_info->pCodecCtx);
        av_close_input_file(st_info->pFormatCtx);
        st_info->pCodecCtx = NULL;
        st_info->pFormatCtx = NULL;
        st_info->pCodec = NULL;
        st_info->width = -1;
        st_info->height = -1;
    }
    return size;
}

__declspec(dllexport) int GetNumberStreams(const char *file) {
    AVFormatContext *pFormatCtx;
    av_register_all();
    // Open video file
    if (av_open_input_file(&pFormatCtx, file, NULL, 0, NULL))
        return -1;  // Couldn't open file

    // Retrieve stream information
    if (av_find_stream_info(pFormatCtx) < 0)
        return -1;  // Couldn't find stream information
    int result = pFormatCtx->nb_streams;
    av_close_input_file(pFormatCtx);
    return result;
}

__declspec(dllexport) long GetNumberVideoFrames(const char *file) {
    long nb_frames = 0L;
    AVFormatContext *pFormatCtx;
    av_register_all();
    // Open video file
    if (av_open_input_file(&pFormatCtx, file, NULL, 0, NULL))
        return -1;  // Couldn't open file

    // Retrieve stream information
    if (av_find_stream_info(pFormatCtx) < 0)
        return -1;  // Couldn't find stream information

    // Find the first video stream
    int videoStream = -1;
    for (unsigned int i = 0; i < pFormatCtx->nb_streams; i++) {
        if (pFormatCtx->streams[i]->codec->codec_type == CODEC_TYPE_VIDEO) {
            videoStream = i;
            break;
        }
    }
    if (videoStream == -1) return -1;  // Didn't find a video stream
    AVStream *str = pFormatCtx->streams[videoStream];
    nb_frames = (long)(pFormatCtx->streams[videoStream]->nb_frames);
    if (nb_frames <= 0) {
        nb_frames =
            (long)av_index_search_timestamp(
                str, str->duration, AVSEEK_FLAG_ANY | AVSEEK_FLAG_BACKWARD) +
            1;
    }
    // Close the video file
    av_close_input_file(pFormatCtx);

    return nb_frames;
}

__declspec(dllexport) float fps(const char *filename) {
    float result = 0;
    AVFormatContext *pFormatCtx;

    // Open video file
    if (av_open_input_file(&pFormatCtx, filename, NULL, 0, NULL))
        return -1;  // Couldn't open file

    // Retrieve stream information
    if (av_find_stream_info(pFormatCtx) < 0)
        return -1;  // Couldn't find stream information

    // Find the first video stream
    int videoStream = -1;
    for (unsigned int i = 0; i < pFormatCtx->nb_streams; i++) {
        if (pFormatCtx->streams[i]->codec->codec_type == CODEC_TYPE_VIDEO) {
            videoStream = i;
            break;
        }
    }
    if (videoStream == -1) return -1;  // Didn't find a video stream

    int num = (pFormatCtx->streams[videoStream]->r_frame_rate).num;
    int den = (pFormatCtx->streams[videoStream]->r_frame_rate).den;
    result = (float)num / (float)den;

    av_close_input_file(pFormatCtx);

    return result;
}

__declspec(dllexport) CImgList<uint8_t> *ph_getKeyFramesFromVideo(
    const char *filename) {
    long N = GetNumberVideoFrames(filename);
    if (N < 0) {
        return NULL;
    }

    float frames_per_sec = (float)(0.5 * fps(filename));
    if (frames_per_sec < 0) {
        return NULL;
    }

    int step = std::round(frames_per_sec);
    long nbframes = (long)(N / step);

    float *dist = (float *)av_malloc((nbframes) * sizeof(float));
    if (!dist) {
        return NULL;
    }
    CImg<float> prev(64, 1, 1, 1, 0);

    VFInfo st_info;
    st_info.filename = strdup(filename);
    st_info.nb_retrieval = 100;
    st_info.step = step;
    st_info.pixelformat = 0;
    st_info.pFormatCtx = NULL;
    st_info.width = -1;
    st_info.height = -1;

    CImgList<uint8_t> *pframelist = new CImgList<uint8_t>();
    if (!pframelist) {
        return NULL;
    }
    int nbread = 0;
    int k = 0;
    do {
        nbread = NextFrames(&st_info, pframelist);
        if (nbread < 0) {
            delete pframelist;
            free(dist);
            return NULL;
        }
        unsigned int i = 0;
        while ((i < (unsigned int)(pframelist->size())) && (k < nbframes)) {
            CImg<uint8_t> current = pframelist->at(i++);
            CImg<float> hist = current.get_histogram(64, 0, 255);
            float d = 0.0;
            dist[k] = 0.0;
            cimg_forX(hist, X) {
                d = hist(X) - prev(X);
                d = (d >= 0) ? d : -d;
                dist[k] += d;
                prev(X) = hist(X);
            }
            k++;
        }
        pframelist->clear();
    } while ((nbread >= st_info.nb_retrieval) && (k < nbframes));
    vfinfo_close(&st_info);

    int S = 10;
    int L = 50;
    int alpha1 = 3;
    int alpha2 = 2;
    int s_begin, s_end;
    int l_begin, l_end;
    uint8_t *bnds = (uint8_t *)av_malloc(nbframes * sizeof(uint8_t));
    if (!bnds) {
        delete pframelist;
        free(dist);
        return NULL;
    }

    int nbboundaries = 0;
    k = 1;
    bnds[0] = 1;
    do {
        s_begin = (k - S >= 0) ? k - S : 0;
        s_end = (k + S < nbframes) ? k + S : nbframes - 1;
        l_begin = (k - L >= 0) ? k - L : 0;
        l_end = (k + L < nbframes) ? k + L : nbframes - 1;

        /* get global average */
        float ave_global, sum_global = 0.0, dev_global = 0.0;
        for (int i = l_begin; i <= l_end; i++) {
            sum_global += dist[i];
        }
        ave_global = sum_global / ((float)(l_end - l_begin + 1));

        /*get global deviation */
        for (int i = l_begin; i <= l_end; i++) {
            float dev = ave_global - dist[i];
            dev = (dev >= 0) ? dev : -1 * dev;
            dev_global += dev;
        }
        dev_global = dev_global / ((float)(l_end - l_begin + 1));

        /* global threshold */
        float T_global = ave_global + alpha1 * dev_global;

        /* get local maximum */
        int localmaxpos = s_begin;
        for (int i = s_begin; i <= s_end; i++) {
            if (dist[i] > dist[localmaxpos]) localmaxpos = i;
        }
        /* get 2nd local maximum */
        int localmaxpos2 = s_begin;
        float localmax2 = 0;
        for (int i = s_begin; i <= s_end; i++) {
            if (i == localmaxpos) continue;
            if (dist[i] > localmax2) {
                localmaxpos2 = i;
                localmax2 = dist[i];
            }
        }
        float T_local = alpha2 * dist[localmaxpos2];
        float Thresh = (T_global >= T_local) ? T_global : T_local;

        if ((dist[k] == dist[localmaxpos]) && (dist[k] > Thresh)) {
            bnds[k] = 1;
            nbboundaries++;
        } else {
            bnds[k] = 0;
        }
        k++;
    } while (k < nbframes - 1);
    bnds[nbframes - 1] = 1;
    nbboundaries += 2;

    int start = 0;
    int end = 0;
    int nbselectedframes = 0;
    do {
        /* find next boundary */
        do {
            end++;
        } while ((bnds[end] != 1) && (end < nbframes));

        /* find min disparity within bounds */
        int minpos = start + 1;
        for (int i = start + 1; i < end; i++) {
            if (dist[i] < dist[minpos]) minpos = i;
        }
        bnds[minpos] = 2;
        nbselectedframes++;
        start = end;
    } while (start < nbframes - 1);

    st_info.nb_retrieval = 1;
    st_info.width = 32;
    st_info.height = 32;
    k = 0;
    do {
        if (bnds[k] == 2) {
            if (ReadFrames(&st_info, pframelist, k * st_info.step,
                           k * st_info.step + 1) < 0) {
                delete pframelist;
                free(dist);
                return NULL;
            }
        }
        k++;
    } while (k < nbframes);
    vfinfo_close(&st_info);

    av_free(bnds);
    bnds = NULL;
    av_free(dist);
    dist = NULL;

    return pframelist;
}
