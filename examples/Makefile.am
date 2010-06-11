INCLUDES = -I$(top_srcdir)/src
noinst_PROGRAMS = buildmvptreedct addmvptreedct querymvptreedct test_texthash test_texthash2

buildmvptreedct_SOURCES = buildmvptree_dctimage.cpp
buildmvptreedct_LDADD = $(top_srcdir)/src/libpHash.la

addmvptreedct_SOURCES = add_mvptree_dct.cpp
addmvptreedct_LDADD = $(top_srcdir)/src/libpHash.la

querymvptreedct_SOURCES = query_mvptree_dct.cpp
querymvptreedct_LDADD = $(top_srcdir)/src/libpHash.la

test_texthash_SOURCES = test_texthash.cpp
test_texthash_LDADD = $(top_srcdir)/src/libpHash.la

test_texthash2_SOURCES = test_texthash2.cpp
test_texthash2_LDADD = $(top_srcdir)/src/libpHash.la

if HAVE_AUDIO_HASH
noinst_PROGRAMS += test_audio build_mvptree_audio add_mvptree_audio query_mvptree_audio

test_audio_SOURCES = test_audiophash.cpp
test_audio_LDADD = $(top_srcdir)/src/libpHash.la

build_mvptree_audio_SOURCES = build_mvptree_audio.cpp
build_mvptree_audio_LDADD = $(top_srcdir)/src/libpHash.la

add_mvptree_audio_SOURCES = add_mvptree_audio.cpp
add_mvptree_audio_LDADD = $(top_srcdir)/src/libpHash.la

query_mvptree_audio_SOURCES = query_mvptree_audio.cpp
query_mvptree_audio_LDADD = $(top_srcdir)/src/libpHash.la
endif

if HAVE_IMAGE_HASH
noinst_PROGRAMS += test_image test_mhimagehash
test_image_SOURCES = test_imagephash.cpp
test_image_LDADD = $(top_srcdir)/src/libpHash.la
test_mhimagehash_SOURCES = test_mhimagehash.cpp
test_mhimagehash_LDADD = $(top_srcdir)/src/libpHash.la
endif

if HAVE_VIDEO_HASH
noinst_PROGRAMS += test_video
test_video_SOURCES = test_dctvideohash.cpp
test_video_LDADD = $(top_srcdir)/src/libpHash.la
endif

