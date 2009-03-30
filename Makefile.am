lib_LTLIBRARIES = libpHash.la
libpHash_la_SOURCES = pHash.cpp audiophash.cpp cimgffmpeg.h

include_HEADERS = pHash.h audiophash.h

bin_PROGRAMS = test_audio test_image test_video
test_audio_SOURCES = test_audiophash.cpp
test_image_SOURCES = test_imagephash.cpp
test_video_SOURCES = test_videophash.cpp

test_audio_LDADD = libpHash.la
test_image_LDADD = libpHash.la
test_video_LDADD = libpHash.la

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = pHash.pc
