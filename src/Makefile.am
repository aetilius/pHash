lib_LTLIBRARIES = libpHash.la
libpHash_la_SOURCES = pHash.cpp cimgffmpeg.cpp ph_fft.c
libpHash_la_LDFLAGS = -no-undefined
include_HEADERS = pHash.h cimgffmpeg.h

if HAVE_AUDIO_HASH
libpHash_la_SOURCES += audiophash.cpp
include_HEADERS += audiophash.h
endif

