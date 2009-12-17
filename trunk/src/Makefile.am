lib_LTLIBRARIES = libpHash.la
libpHash_la_SOURCES = pHash.cpp ph_fft.c
libpHash_la_LDFLAGS = -no-undefined
include_HEADERS = pHash.h cimgffmpeg.h audiophash.h ph_fft.h

if HAVE_AUDIO_HASH
libpHash_la_SOURCES += audiophash.cpp audiophash.h
endif

