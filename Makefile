# pHash project makefile v0.0.1

CC      = g++
CCFLAGS = -Wall -ffast-math -O2
OUTFILE = pHash
TESTFILE = test_main.cpp
TEST2FILE = dct_image_main.cpp 
LIBS = -lm -lpthread -ljpeg -lpHash
X11LIBS = -lX11 -lXext -lXrandr 
FFMPEGLIBS = -lavformat -lavcodec -lavutil -lswscale
X11LIBDIRS = -L/usr/X11R6/lib64
FFMPEGLIBDIRS = -L/usr/local/lib
LIBDIRS = -L.
FFMPEGINCLUDEDIRS = -I/usr/include/ffmpeg
CIMGDEFINES = -Dcimg_use_jpeg -Dcimg_display=0 -Dcimg_debug=0
X11DEFINES = -Dcimg_use_xshm -Dcimg_use_xrandr



test_rash_image : pHash.so
	$(CC) $(CCFLAGS) $(TESTFILE) $(CIMGDEFINES) -o$(OUTFILE) $(LIBDIRS) $(LIBS)

test_dct_image: pHash.so
	$(CC) $(CCFLAGS) $(TEST2FILE) $(CIMGDEFINES) -opHash2  $(LIBDIRS) $(LIBS) 

test_dct_video: pHash.so
	$(CC) $(CCFLAGS) dct_video_main.cpp $(CIMGDEFINES) $(FFMPEGINCLUDEDIRS) -opHash3 $(LIBDIRS) $(LIBS) $(FFMPEGLIBDIRS) $(FFMPEGLIBS)

pHash.a : pHash.o
	ar rcs libpHash.a pHash.o

pHash.so : pHash.o
	$(CC) -shared $(CCFLAGS) $(CIMGDEFINES) pHash.o -Wl,-soname -Wl,libpHash.so.0.2 -olibpHash.so.0.2 $(FFMPEGLIBS) 
	ln -sf libpHash.so.0.2 libpHash.so

pHash.o : pHash.cpp pHash.h
	$(CC) -fPIC -DPIC $(CIMGDEFINES) $(CCFLAGS) -c pHash.cpp

clean :
	rm pHash2 pHash.so.* 
