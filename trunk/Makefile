PHASH_VERSION = 004
PREFIX = /usr/local
CC      = g++
CCFLAGS = -Wall -O2 -ffast-math
OUTFILE = pHash
TESTFILE = test_main.cpp
TEST2FILE = dct_image_main.cpp 
LIBS = -lm -lpthread -ljpeg -lpng -lpHash
FFTW3LIBS = -lfftw3
FFMPEGLIBS = -lavformat -lavcodec -lavutil -lswscale
FFMPEGLIBDIRS = -L/usr/local/lib
LIBDIRS = -L.
FFMPEGINCLUDEDIRS = -I/usr/include/ffmpeg -I/usr/local/include/
CIMGDEFINES = -Dcimg_use_jpeg -Dcimg_use_png -Dcimg_display=0 -Dcimg_debug=0 -DPHASH_VERSION=$(PHASH_VERSION)


install : pHash.so
	cp libpHash.so.0.4 $(PREFIX)/lib
	cp *.h $(PREFIX)/include
	ln -sf $(PREFIX)/lib/libpHash.so.0.4 $(PREFIX)/lib/libpHash.so
	ldconfig -n $(PREFIX)/lib

test_audio: pHash.so
	$(CC) $(CCFLAGS) test_audiophash.cpp $(FFMPEGINCLUDEDIRS) -otest_audiophash $(LIBDIRS) $(LIBS) $(FFMPEGLIBDIRS) $(FFMPEGLIBS) $(FFTW3LIBS)

test_image: pHash.so
	$(CC) $(CCFLAGS) test_imagephash.cpp $(FFMPEGINCLUDEDIRS) -otest_imagephash $(LIBDIRS) $(LIBS) $(FFMPEGLIBDIRS) $(FFMPEGLIBS) $(FFTW3LIBS)

test_video: pHash.so
	$(CC) $(CCFLAGS) test_videophash.cpp $(FFMPEGINCLUDEDIRS) -otest_videophash $(LIBDIRS) $(LIBS) $(FFMPEGLIBDIRS) $(FFMPEGLIBS) $(FFTW3LIBS)

pHash.a : pHash.o audiophash.o
	ar rcs libpHash.a *.o

pHash.so : pHash.o audiophash.o
	$(CC) -shared $(CCFLAGS) $(CIMGDEFINES) *.o  -Wl,-soname -Wl,libpHash.so.0.4 -olibpHash.so.0.4 $(FFMPEGLIBDIRS) $(FFMPEGLIBS) $(FFTW3LIBS)
	ln -sf libpHash.so.0.4 libpHash.so

pHash.o : pHash.cpp pHash.h 
	$(CC) -fPIC -DPIC $(CIMGDEFINES) $(CCFLAGS) $(FFMPEGINCLUDEDIRS) -c pHash.cpp


audiophash.o: audiophash.cpp audiophash.h
	$(CC) -fPIC -DPIC $(CCFLAGS) $(FFMPEGINCLUDEDIRS)  -c audiophash.cpp

clean :
	rm -f pHash pHash3 pHash2 libpHash.* *.o test_rash_image audiophash
