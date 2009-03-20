#include <jni.h>

#ifndef _Included_pHash
#define _Included_pHash
#ifdef __cplusplus
extern "C" {
#endif

JNIEXPORT jlong JNICALL Java_pHash_imageHash
  (JNIEnv *, jclass, jstring);

JNIEXPORT jlong JNICALL Java_pHash_videoHash
  (JNIEnv *, jclass, jstring);

JNIEXPORT jintArray JNICALL Java_pHash_audioHash
  (JNIEnv *, jclass, jstring);

JNIEXPORT jint JNICALL Java_pHash_distance
  (JNIEnv *, jclass, jlong, jlong);

#ifdef __cplusplus
}
#endif
#endif
