lib_LTLIBRARIES = libpHash.la
libpHash_la_SOURCES = pHash.cpp
libpHash_la_LDFLAGS = -no-undefined
include_HEADERS = pHash.h
bin_PROGRAMS =

EXAMPLES = 
if HAVE_AUDIO_HASH
bin_PROGRAMS += test_audio

libpHash_la_SOURCES += audiophash.cpp
test_audio_SOURCES = test_audiophash.cpp
test_audio_LDADD = libpHash.la
endif

if HAVE_IMAGE_HASH
bin_PROGRAMS += test_image
test_image_SOURCES = test_imagephash.cpp
test_image_LDADD = libpHash.la
endif

if HAVE_VIDEO_HASH
bin_PROGRAMS += test_video
test_video_SOURCES = test_videophash.cpp
test_video_LDADD = libpHash.la

endif

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = pHash.pc
