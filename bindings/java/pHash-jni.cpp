
#include "pHash-jni.h"
#include <jni.h>
#include "pHash.h"
#include "audiophash.h"
//jmethodid bitSet_set1;
//jclass bitSet_class1;
//jmethodid bitSet_constructor1;

JNIEXPORT void JNICALL Java_pHashd_initIDs
  (JNIEnv *e, jclass cl)
{
//    bitSet_class1 = e->FindClass("java/util/BitSet");
 //   bitSet_constructor1 = e->GetMethodID(bitSet_class1, "<init>", "(I)V");

   // bitSet_set1 = e->GetMethodID(bitSet_class1, "set", "(I)V");


}
JNIEXPORT jobject JNICALL Java_pHash_imageHash
  (JNIEnv *e, jclass cl, jstring f)
{
    
    static jclass bitSet_class = NULL;
    static jmethodID bitSet_constructor = NULL;
    static jmethodID bitSet_set = NULL;

    if(bitSet_class == NULL)
    {
        bitSet_class = e->FindClass("java/util/BitSet");
        if(bitSet_class == NULL)
            return NULL;
    }

    if(bitSet_constructor == NULL)
    {
        bitSet_constructor = e->GetMethodID(bitSet_class, "<init>", "(I)V");
        if(bitSet_constructor == NULL)
            return NULL;
    }
    if(bitSet_set == NULL)
    {
        bitSet_set = e->GetMethodID(bitSet_class, "set", "(I)V");
        if(bitSet_set == NULL)
            return NULL;
    }

	const char *file = e->GetStringUTFChars(f,0);
    jobject ret = NULL;

	ret = e->NewObject(bitSet_class, bitSet_constructor, 64);

    ulong64 hash;
    ph_dct_imagehash(file, hash);
    jint i = 0;
    while(hash)
    {
        e->CallVoidMethod(bitSet_class, bitSet_set, i);
        ++i;
        hash &= hash - 1;
    }
	e->ReleaseStringUTFChars(f,file);
	
	return ret;
	
}

JNIEXPORT jintArray JNICALL Java_pHash_audioHash
  (JNIEnv *e, jclass cl, jstring f)
{
    
    static jclass bitSet_class = NULL;
    static jmethodID bitSet_constructor = NULL;
    static jmethodID bitSet_set = NULL;

    if(bitSet_class == NULL)
    {
        bitSet_class = e->FindClass("java/util/BitSet");
        if(bitSet_class == NULL)
            return NULL;
    }

    if(bitSet_constructor == NULL)
    {
        bitSet_constructor = e->GetMethodID(bitSet_class, "<init>", "(I)V");
        if(bitSet_constructor == NULL)
            return NULL;
    }
    if(bitSet_set == NULL)
    {
        bitSet_set = e->GetMethodID(bitSet_class, "set", "(I)V");
        if(bitSet_set == NULL)
            return NULL;
    }

	const char *file = e->GetStringUTFChars(f,0);
    jintArray ret = NULL;
    const float threshold = 0.30;        //ber threshold (0.25-0.35)
    const int block_size = 256;          //number of frames to compare at a time
    const int sr = 8000;                 //sample rate to convert the stream
    const int channels = 1;              //number of channels to convert stream
    int N;
    int nbframes;
    float *buf = NULL;
    unsigned int *hash = NULL;

    buf = ph_readaudio(file,sr,channels,N); 
    e->ReleaseStringUTFChars(f,file);
    hash = ph_audiohash(buf,N,sr,nbframes);
    free(buf);
    ret = e->NewIntArray(nbframes);

    e->SetIntArrayRegion(ret, 0, nbframes, hash);

    free(hash);

	
	return ret;
	
}