
#include "pHash-jni.h"
#include <jni.h>
#include "pHash.h"
#include "audiophash.h"

JNIEXPORT jint JNICALL Java_pHash_distance
  (JNIEnv *e, jclass cl, jlong hash1, jlong hash2)
{
    
	return ph_hamming_distance(hash1, hash2);
	
}

JNIEXPORT jlong JNICALL Java_pHash_imageHash
  (JNIEnv *e, jclass cl, jstring f)
{
    
	const char *file = e->GetStringUTFChars(f,0);
    
    ulong64 hash;
    ph_dct_imagehash(file, hash);
    e->ReleaseStringUTFChars(f,file);

	return hash;
	
}

JNIEXPORT jlong JNICALL Java_pHash_videoHash
  (JNIEnv *e, jclass cl, jstring f)
{
    
	const char *file = e->GetStringUTFChars(f,0);
    
    ulong64 hash;
    ph_dct_videohash(file, hash);
    e->ReleaseStringUTFChars(f,file);

	return hash;
	
}
JNIEXPORT jintArray JNICALL Java_pHash_audioHash
  (JNIEnv *e, jclass cl, jstring f)
{
    
    jintArray ret = NULL;
    const int sr = 8000;
    const int channels = 1;
    int N;
    int nbframes;
    float *buf = NULL;
    unsigned int *hash = NULL;
	const char *file = e->GetStringUTFChars(f,0);
    buf = ph_readaudio(file,sr,channels,N); 
    e->ReleaseStringUTFChars(f,file);
    hash = ph_audiohash(buf,N,sr,nbframes);
    free(buf);
    ret = e->NewIntArray(nbframes);

    e->SetIntArrayRegion(ret, 0, nbframes, hash);

    free(hash);

	return ret;
	
}