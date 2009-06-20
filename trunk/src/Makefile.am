lib_LTLIBRARIES = libpHash.la
libpHash_la_SOURCES = pHash.cpp
libpHash_la_LDFLAGS = -no-undefined
include_HEADERS = pHash.h

if HAVE_AUDIO_HASH
libpHash_la_SOURCES += audiophash.cpp
endif

