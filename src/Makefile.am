lib_LTLIBRARIES = libpHash.la
libpHash_la_SOURCES = pHash.cpp
libpHash_la_LDFLAGS = -no-undefined
include_HEADERS = pHash.h

if HAVE_AUDIO_HASH
libpHash_la_SOURCES += audiophash.cpp ph_fft.c
include_HEADERS += audiophash.h ph_fft.h
endif

if HAVE_VIDEO_HASH
libpHash_la_SOURCES += cimgffmpeg.cpp cimgffmpeg.h
endif

