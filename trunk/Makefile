# pHash project makefile v0.0.1

CC      = g++
CCFLAGS = -Wall -ffast-math -O2
OUTFILE = pHash
TESTFILE = test_main.cpp
TEST2FILE = dct_image_main.cpp 
LIBS = -lm -lpthread -ljpeg -lpHash
X11LIBS = -lX11 -lXext -lXrandr 
X11LIBDIRS = -L/usr/X11R6/lib64
LIBDIRS = -L.
CIMGDEFINES = -Dcimg_use_jpeg -Dcimg_display=0 -Dcimg_debug=0
X11DEFINES = -Dcimg_use_xshm -Dcimg_use_xrandr

test : pHash.so
	$(CC) $(CCFLAGS) $(TESTFILE) $(CIMGDEFINES) -o$(OUTFILE) $(LIBDIRS) $(LIBS)

test2: pHash.so
	$(CC) $(CCFLAGS) $(TEST2FILE) $(CIMGDEFINES) -opHash2  $(LIBDIRS) $(LIBS) 

pHash.so : pHash.o
	$(CC) -shared $(CCFLAGS) $(CIMGDEFINES) pHash.o -Wl,-soname -Wl,libpHash.so.0.2 -olibpHash.so.0.2  
	ln -sf libpHash.so.0.2 libpHash.so

pHash.o : pHash.cpp pHash.h
	$(CC) -fPIC -DPIC $(CIMGDEFINES) $(CCFLAGS) -c pHash.cpp

clean :
	rm pHash2 pHash.so.* 
