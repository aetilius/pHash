
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

	if(hash1 == NULL || hash2 == NULL)
	{
		return (jdouble)-1.0;
	}
	const float threshold = 0.30;
    	const int block_size = 256;
    	int Nc;
	jdouble maxC = 0.0;
	double *pC; 
	jint hash1_len, hash2_len;
	hash1_len = e->GetArrayLength(hash1);
	hash2_len = e->GetArrayLength(hash2);
	if(hash1_len <= 0 || hash2_len <= 0)
	{
		return (jdouble)-1.0;
	}
	uint32_t *hash1_n, *hash2_n;
	
	hash1_n = (uint32_t*)e->GetIntArrayElements(hash1, 0);
	hash2_n = (uint32_t*)e->GetIntArrayElements(hash2, 0);
	pC = ph_audio_distance_ber(hash1_n, hash1_len, hash2_n, hash2_len, threshold, block_size, Nc);
	e->ReleaseIntArrayElements(hash1, (jint*)hash1_n, 0);
	e->ReleaseIntArrayElements(hash2, (jint*)hash2_n, 0);
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
	if(!buf) 
	{
    		e->ReleaseStringUTFChars(f,file);
		return NULL;
	}
    e->ReleaseStringUTFChars(f,file);
    hash = ph_audiohash(buf,N,sr,nbframes);
	if(!hash) 
	{
    		free(buf);
		return NULL;
	}
	free(buf);
    ret = e->NewIntArray(nbframes);

    e->SetIntArrayRegion(ret, 0, nbframes, (jint *)hash);

    free(hash);

	return ret;
	
}
