INCLUDES = -I$(top_srcdir)/src
noinst_PROGRAMS = test_texthash test_texthash2

test_texthash_SOURCES = test_texthash.cpp
test_texthash_LDADD = $(top_builddir)/src/libpHash.la

test_texthash2_SOURCES = test_texthash2.cpp
test_texthash2_LDADD = $(top_builddir)/src/libpHash.la

if HAVE_AUDIO_HASH
noinst_PROGRAMS += test_audio 

test_audio_SOURCES = test_audiophash.cpp
test_audio_LDADD = $(top_builddir)/src/libpHash.la

endif

if HAVE_IMAGE_HASH
noinst_PROGRAMS += test_image test_mhimagehash

test_image_SOURCES = test_imagephash.cpp
test_image_LDADD = $(top_builddir)/src/libpHash.la
test_mhimagehash_SOURCES = test_mhimagehash.cpp
test_mhimagehash_LDADD = $(top_builddir)/src/libpHash.la
endif

if HAVE_VIDEO_HASH
noinst_PROGRAMS += test_video
test_video_SOURCES = test_dctvideohash.cpp
test_video_LDADD = $(top_builddir)/src/libpHash.la
endif
