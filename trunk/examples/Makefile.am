INCLUDES = -I$(top_srcdir)/src
bin_PROGRAMS = build_mvptree add_mvptree query_mvptree

build_mvptree_SOURCES = build_mvptree.cpp
build_mvptree_LDADD = $(top_srcdir)/src/libpHash.la

add_mvptree_SOURCES = add_mvptree.cpp
add_mvptree_LDADD = $(top_srcdir)/src/libpHash.la

query_mvptree_SOURCES = query_mvptree.cpp
query_mvptree_LDADD = $(top_srcdir)/src/libpHash.la

if HAVE_AUDIO_HASH
bin_PROGRAMS += test_audio build_mvptree_audio add_mvptree_audio query_mvptree_audio

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
bin_PROGRAMS += test_image
test_image_SOURCES = test_imagephash.cpp
test_image_LDADD = $(top_srcdir)/src/libpHash.la
endif

if HAVE_VIDEO_HASH
bin_PROGRAMS += test_video
test_video_SOURCES = test_videophash.cpp
test_video_LDADD = $(top_srcdir)/src/libpHash.la

endif

