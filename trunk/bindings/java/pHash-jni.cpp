
#include "pHash-jni.h"
#include <jni.h>
#include "pHash.h"
#include "audiophash.h"

JNIEXPORT jint JNICALL Java_pHash_imageDistance
  (JNIEnv *e, jclass cl, jlong hash1, jlong hash2)
{
    
	return ph_hamming_distance((ulong64)hash1, (ulong64)hash2);
	
}
JNIEXPORT jdouble JNICALL Java_pHash_audioDistance
  (JNIEnv *e, jclass cl, jintArray hash1, jintArray hash2)
{

	const float threshold = 0.30;
    	const int block_size = 256;
    	int Nc;
	jdouble maxC = 0.0;
	double *pC; 
	jint hash1_len, hash2_len;
	hash1_len = e->GetArrayLength(hash1);
	hash2_len = e->GetArrayLength(hash2);
	if(hash1_len <= 0 || hash2_len <= 0)
		return -1;
	jint *hash1_n, *hash2_n;
	hash1_n = (jint *)malloc(sizeof(jint)*hash1_len);
	hash2_n = (jint *)malloc(sizeof(jint)*hash2_len);
	
	e->GetIntArrayRegion(hash1, 0, hash1_len, hash1_n);
	e->GetIntArrayRegion(hash2, 0, hash2_len, hash2_n);
	pC = ph_audio_distance_ber((uint32_t*)hash1_n, hash1_len, (uint32_t*)hash2_n, hash2_len, threshold, block_size, Nc);
	free(hash1_n);
	free(hash2_n);
	maxC = 0.0;
        for (int j=0;j<Nc;j++){
            if (pC[j] > maxC){
                maxC = pC[j];
            }
	}
	delete[] pC;
	return maxC;
	
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

    e->SetIntArrayRegion(ret, 0, nbframes, (jint *)hash);

    free(hash);

	return ret;
	
}
