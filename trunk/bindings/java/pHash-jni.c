#include "pHash-jni.h"
#include "pHash.h"

JNIEXPORT jobject JNICALL Java_pHash_calculateHash
  (JNIEnv *e, jclass class, jstring f, jint hashtype)
{

	const char *file = (*e)->GetStringUTFChars(e,f,0);

	jclass c = (*e)->FindClass(e,"java/util/BitSet");
	jmethodID mid = (*e)->GetMethodID(e, c,"<init>", "(I)V");
	jclass ret = (*e)->NewObject(e, c, mid, 64);
	
	
	(*e)->ReleaseStringUTFChars(e,f,file);
	(*e)->DeleteLocalRef(e,c);
	
	return ret;
	
}
