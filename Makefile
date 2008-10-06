# pHash project makefile v0.0.1

CC      = g++
CCFLAGS = -Wall
OUTFILE = pHash
TESTFILE = test_main.cpp
LIBS = -lm -lpthread -ljpeg -lpHash
LIBDIRS = -L.
CIMGDEFINES = -Dcimg_use_jpeg -Dcimg_display=0

test : pHash.so
	$(CC) $(CCFLAGS) $(TESTFILE) $(CIMGDEFINES) -o$(OUTFILE) $(LIBDIRS) $(LIBS)

pHash.so : pHash.o
	$(CC) -shared $(CIMGDEFINES) pHash.o -Wl,-soname -Wl,libpHash.so.0.2 -olibpHash.so.0.2  
	ln -sf libpHash.so.0.2 libpHash.so

pHash.o : pHash.cpp pHash.h
	$(CC) -fPIC -DPIC $(CIMGDEFINES) $(CCFLAGS) -c pHash.cpp

clean :
	rm pHash pHash.o
