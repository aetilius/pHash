lib_LTLIBRARIES = libpHash-jni.la
libpHash_jni_la_SOURCES = pHash-jni.cpp
include_HEADERS = pHash-jni.h

clean:
	rm -f *.class *.so *.la *.lo

all: pHash.class
pHash.java: libpHash_jni.la

pHash.class: pHash.java
	$(GCJ) -C pHash.java 

libpHash_jni_la_LIBADD = ../../libpHash.la
