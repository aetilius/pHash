#include "pHash-jni.h"
#include "pHash.h"

JNIEXPORT jobject JNICALL Java_pHash_calculateHash
  (JNIEnv *e, jclass cl, jstring f, jint hashtype)
{

	const char *file = e->GetStringUTFChars(f,0);

	jclass c = e->FindClass("java/util/BitSet");
	jmethodID mid = e->GetMethodID(c,"<init>", "(I)V");
	jobject ret = e->NewObject(c, mid, 64);
		
	
	e->ReleaseStringUTFChars(f,file);
	e->DeleteLocalRef(c);
	
	return ret;
	
}
