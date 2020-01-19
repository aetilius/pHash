#include "pHash-jni.h"
#include <pHash.h>
#include <audiophash.h>


static jint JNI_VERSION = JNI_VERSION_1_8;
	
// classes
static jclass bmbHashClass;
static jclass digestClass;
static jclass featuresClass;
static jclass mhHashClass;
static jclass txtHashPointClass;
static jclass txtMatchClass;

// ctor ids
static jmethodID bmbHashCtor;
static jmethodID digestCtor;
static jmethodID featuresCtor;
static jmethodID mhHashCtor;
static jmethodID txtHashPointCtor;
static jmethodID txtMatchCtor;

// image getters
static jmethodID bmbHash_getHash_bytearray;
static jmethodID digest_getCoeffs_bytearray;
static jmethodID features_getFeatures_doublearray;
static jmethodID mhHash_getHash_bytearray;
	
// text getters
static jmethodID txtHashPoint_getHash_long;
static jmethodID txtHashPoint_getIndex_long;
static jmethodID txtMatch_getFirstIndex_long;
static jmethodID txtMatch_getSecondIndex_long;
static jmethodID txtMatch_getLength_long;

// text setters
static jmethodID txtHashPoint_setDataMembers_longlong;
static jmethodID txtMatch_setDataMembers_longlonglong;
	
	
jint JNI_OnLoad(JavaVM *vm, void *reserved){
	JNIEnv *env = NULL;
	if (vm->GetEnv(reinterpret_cast<void**>(&env), JNI_VERSION) != JNI_OK)
		return JNI_ERR;

	// global references to classes
	jclass tmpBMBHashClass = env->FindClass("org/phash/ni/image/BMBHash");
	bmbHashClass = (jclass)env->NewGlobalRef(tmpBMBHashClass);
	env->DeleteLocalRef(tmpBMBHashClass);
	
	jclass tmpDigestClass = env->FindClass("org/phash/ni/image/Digest");
	digestClass = (jclass)env->NewGlobalRef(tmpDigestClass);
	env->DeleteLocalRef(tmpDigestClass);
	
	jclass tmpFeaturesClass = env->FindClass("org/phash/ni/image/Features");
	featuresClass = (jclass)env->NewGlobalRef(tmpFeaturesClass);
	env->DeleteLocalRef(tmpFeaturesClass);

	jclass tmpmhHashClass = env->FindClass("org/phash/ni/image/MHash");
	mhHashClass = (jclass)env->NewGlobalRef(tmpmhHashClass);
	env->DeleteLocalRef(mhHashClass);

	jclass tmpTxtHashPointClass = env->FindClass("org/phash/ni/text/TxtHashPoint");
	txtHashPointClass = (jclass)env->NewGlobalRef(tmpTxtHashPointClass);
	env->DeleteLocalRef(tmpTxtHashPointClass);

	jclass tmpTxtMatchClass = env->FindClass("org/phash/ni/text/TxtMatchPoint");
	txtMatchClass = (jclass)env->NewGlobalRef(tmpTxtMatchClass);
	env->DeleteLocalRef(tmpTxtMatchClass);

	// ctor methodIDs 
	bmbHashCtor      = env->GetMethodID(bmbHashClass, "<init>", "([B)V");
	digestCtor       = env->GetMethodID(digestClass, "<init>", "([B)V");
	featuresCtor     = env->GetMethodID(featuresClass, "<init>", "([D)V");
	mhHashCtor       = env->GetMethodID(mhHashClass, "<init>", "([B)V");
	txtHashPointCtor = env->GetMethodID(txtHashPointClass, "<init>", "(JJ)V");
	txtMatchCtor     = env->GetMethodID(txtMatchClass, "<init>", "(JJJ)V");

	// getter image methodIDs
	bmbHash_getHash_bytearray        = env->GetMethodID(bmbHashClass, "getHash", "()[B");
	digest_getCoeffs_bytearray       = env->GetMethodID(digestClass, "getCoeffs", "()[B");
	features_getFeatures_doublearray = env->GetMethodID(featuresClass, "getFeatures", "()[D");
	mhHash_getHash_bytearray         = env->GetMethodID(mhHashClass, "getHash", "()[B");
	
	// getter text methodIDs
	txtHashPoint_getHash_long    = env->GetMethodID(txtHashPointClass, "getHash", "()J");
	txtHashPoint_getIndex_long   = env->GetMethodID(txtHashPointClass, "getIndex", "()J");
	txtMatch_getFirstIndex_long  = env->GetMethodID(txtMatchClass, "getFirstIndex", "()J");
	txtMatch_getSecondIndex_long = env->GetMethodID(txtMatchClass, "getSecondIndex", "()J");
	txtMatch_getLength_long      = env->GetMethodID(txtMatchClass, "getLength", "()J");

	// setter text methodIDs
	txtHashPoint_setDataMembers_longlong = env->GetMethodID(txtHashPointClass, "setDataMembers", "(JJ)V");
	txtMatch_setDataMembers_longlonglong = env->GetMethodID(txtMatchClass, "setDataMembers", "(JJJ)V");
		
}

void JNI_OnUnload(JavaVM *vm, void *reserved){
	JNIEnv *env = NULL;
	if (vm->GetEnv(reinterpret_cast<void**>(&env), JNI_VERSION) != JNI_OK)
		return;
	
	env->DeleteGlobalRef(bmbHashClass);
	env->DeleteGlobalRef(digestClass);
	env->DeleteGlobalRef(featuresClass);
	env->DeleteGlobalRef(mhHashClass);
	env->DeleteGlobalRef(txtHashPointClass);
	env->DeleteGlobalRef(txtMatchClass);
	return;
}
	
#ifdef USE_IMAGE_HASH

JNIEXPORT jlong JNICALL Java_org_phash_ni_image_ImageHash_dctImageHash
(JNIEnv *env, jobject thisObj, jstring file){

	const char *filename = env->GetStringUTFChars(file, NULL);
	if (filename == NULL)
		return -1;
	
	ulong64 hashval;
	if (ph_dct_imagehash(filename, hashval) < 0){
		env->ReleaseStringUTFChars(file, filename);
		return -1;
	}
	
	env->ReleaseStringUTFChars(file, filename);

	/* cast to jlong while preserving bit pattern */
	jlong result = ((union { ulong64 i; jlong result;}){.i=hashval}).result;
	return result;
}

JNIEXPORT jint JNICALL Java_org_phash_ni_image_ImageHash_dctImageDistance
(JNIEnv *env, jobject thisObj, jlong hash1, jlong hash2){

	/* cast hash args to ulong64 while preserving bit pattern */
	ulong64 h1 = ((union {jlong j; ulong64 h1;}){.j=hash1}).h1;
	ulong64 h2 = ((union {jlong j; ulong64 h2;}){.j=hash2}).h2;
	jint d = (jint)ph_hamming_distance(h1, h2);;
	return d;
}

JNIEXPORT jobject JNICALL Java_org_phash_ni_image_ImageHash_bmbImageHash
(JNIEnv *env, jobject thisObj, jstring file){

	const char *filename = env->GetStringUTFChars(file, NULL);

	BMBHash bmbhash;
	if (ph_bmb_imagehash(filename, bmbhash) < 0){
		env->ReleaseStringUTFChars(file, filename);
		return NULL;
	}

	/* create new java array and BMBHash object */
	jbyteArray arr = env->NewByteArray((jsize)bmbhash.bytelength);
	env->SetByteArrayRegion(arr, 0, (jsize)bmbhash.bytelength, (jbyte*)bmbhash.hash);
	
	jobject bmbHashObj = env->NewObject(bmbHashClass, bmbHashCtor, arr);
	
	env->ReleaseStringUTFChars(file, filename);

	return bmbHashObj;
}

JNIEXPORT jint JNICALL Java_org_phash_ni_image_ImageHash_bmbDistance
(JNIEnv *env, jobject thisObj, jobject hash1, jobject hash2){
	if ((env->IsInstanceOf(hash1, bmbHashClass) == JNI_FALSE) ||
		(env->IsInstanceOf(hash2, bmbHashClass) == JNI_FALSE)){
		return -1;
	}

	jbyteArray data1 = reinterpret_cast<jbyteArray>(env->CallObjectMethod(hash1, bmbHash_getHash_bytearray));
	uint32_t n1 = (uint32_t)env->GetArrayLength(data1);
	
	jbyteArray data2 = reinterpret_cast<jbyteArray>(env->CallObjectMethod(hash2, bmbHash_getHash_bytearray));
	uint32_t n2 = (uint32_t)env->GetArrayLength(data2);

	jbyte *pdata1 = env->GetByteArrayElements(data1, NULL);
	jbyte *pdata2 = env->GetByteArrayElements(data2, NULL);
	
	BMBHash h1;
	h1.hash = reinterpret_cast<uint8_t*>(pdata1);
	h1.bytelength = static_cast<uint32_t>(n1);

	BMBHash h2;
	h2.hash = reinterpret_cast<uint8_t*>(pdata2);
	h2.bytelength = static_cast<uint32_t>(n2);
	
	jint d = (jint)ph_bmb_distance(h1, h2);

	env->ReleaseByteArrayElements(data1, pdata1, 0);
	env->ReleaseByteArrayElements(data2, pdata2, 0);
	
	return d;
}

JNIEXPORT jobject JNICALL Java_org_phash_ni_image_ImageHash_radialFeatureVector
(JNIEnv *env, jobject thisObj, jstring file, jint nAngles, jdouble sigma, jdouble gamma){
	const char *filename = env->GetStringUTFChars(file, NULL);

	Features fv;
	if (ph_feature_vector(filename, (int)nAngles, (double)sigma, (double)gamma, fv) < 0){
		env->ReleaseStringUTFChars(file, filename);
		return NULL;
	}

	jdoubleArray arr = env->NewDoubleArray(fv.n_features);
	env->SetDoubleArrayRegion(arr, 0, (jsize)fv.n_features, (jdouble*)fv.features);
	
	jobject fvObj = env->NewObject(featuresClass, featuresCtor, arr);
	
	env->ReleaseStringUTFChars(file, filename);
	return fvObj;
}


JNIEXPORT jobject JNICALL Java_org_phash_ni_image_ImageHash_imageDigest
(JNIEnv *env, jobject thisObj, jstring file, jdouble sigma, jdouble gamma, jint nAngles){
	const char *filename = env->GetStringUTFChars(file, NULL);

	Digest digest;
	if (ph_image_digest(filename, (double)sigma, (double)gamma, digest, (int)nAngles) < 0){
		env->ReleaseStringUTFChars(file, filename);
		return NULL;
	}

	env->ReleaseStringUTFChars(file, filename);

	jbyteArray arr = env->NewByteArray((jsize)digest.n_coeffs);
	env->SetByteArrayRegion(arr, 0, (jsize)digest.n_coeffs, (jbyte*)digest.coeffs);
	jobject digestObj = env->NewObject(digestClass, digestCtor, arr);

	return digestObj;
}

JNIEXPORT jdouble JNICALL Java_org_phash_ni_image_ImageHash_peakCrossCorr
(JNIEnv *env, jobject thisObj, jobject digest1, jobject digest2){
	if ((env->IsInstanceOf(digest1, digestClass) == JNI_ERR) ||
		(env->IsInstanceOf(digest2, digestClass) == JNI_ERR))
		return -1;

	jbyteArray data1 = reinterpret_cast<jbyteArray>(env->CallObjectMethod(digest1, digest_getCoeffs_bytearray));
	jbyte *pdata1 = env->GetByteArrayElements(data1, NULL);
	
	jsize n1 = env->GetArrayLength(data1);

	jbyteArray data2 = reinterpret_cast<jbyteArray>(env->CallObjectMethod(digest2, digest_getCoeffs_bytearray));
	jbyte *pdata2 = env->GetByteArrayElements(data2, NULL);
	jsize n2 = env->GetArrayLength(data2);

	Digest x;
	x.coeffs = reinterpret_cast<uint8_t*>(pdata1);
	x.n_coeffs = (uint32_t)n1;

	Digest y;
	y.coeffs = reinterpret_cast<uint8_t*>(pdata2);
	y.n_coeffs = (uint32_t)n2;

	env->ReleaseByteArrayElements(data1, pdata1, 0);
	env->ReleaseByteArrayElements(data2, pdata2, 0);

	double pcc;
	if (ph_peakcrosscorr(x, y, pcc) < 0){
		return -1;
	}

	ph_free_digest(x);
	ph_free_digest(y);
	
	return (jdouble)pcc;
}

JNIEXPORT jobject JNICALL Java_org_phash_ni_image_ImageHash_mhImageHash
(JNIEnv *env, jobject thisObj, jstring file, jfloat alpha, jfloat level){
	const char *filename = env->GetStringUTFChars(file, NULL);

	int n;
	uint8_t *mhhash = ph_mh_imagehash(filename, n, (float)alpha, (float)level);
	if (mhhash == NULL){
		env->ReleaseStringUTFChars(file, filename);
		return NULL;
	}

	jbyteArray arr = env->NewByteArray((jsize)n);
	env->SetByteArrayRegion(arr, 0, (jsize)n, (jbyte*)mhhash);

	jobject mhObj = env->NewObject(mhHashClass, mhHashCtor, arr);
	env->ReleaseStringUTFChars(file, filename);
	return mhObj;
}

JNIEXPORT jint JNICALL Java_org_phash_ni_image_ImageHash_mhImageDistance
(JNIEnv *env, jobject thisObj, jobject hash1, jobject hash2){
	if ((env->IsInstanceOf(hash1, mhHashClass) == JNI_ERR) ||
		(env->IsInstanceOf(hash2, mhHashClass) == JNI_ERR)){
		return -1;
	}

	jbyteArray data1 = reinterpret_cast<jbyteArray>(env->CallObjectMethod(mhHashClass, mhHash_getHash_bytearray));
	jbyte *pdata1 = env->GetByteArrayElements(data1, NULL);
	jsize n1 = env->GetArrayLength(data1);

	jbyteArray data2 = reinterpret_cast<jbyteArray>(env->CallObjectMethod(mhHashClass, mhHash_getHash_bytearray));
	jbyte *pdata2 = env->GetByteArrayElements(data2, NULL);
	jsize n2 = env->GetArrayLength(data2);
	
	int d = ph_hammingdistance2((uint8_t*)pdata1, (int)n1, (uint8_t*)pdata2, (int)n2);

	env->ReleaseByteArrayElements(data1, pdata1, 0);
	env->ReleaseByteArrayElements(data2, pdata2, 0);
	
	return (jint)d;

}

#endif /* USE_IMAGE_HASH */




	
#ifdef USE_AUDIO_HASH

JNIEXPORT jint JNICALL Java_org_phash_ni_audio_AudioHash_countAudioSamples
(JNIEnv *env, jobject thisObj, jstring file, jint sr, jint nChannels){
	const char *filename = env->GetStringUTFChars(file, NULL);
	int n_samples = ph_count_samples(filename, (int)sr, (int)nChannels);
	env->ReleaseStringUTFChars(file, filename);
	return (jint)n_samples;
}

JNIEXPORT jfloatArray JNICALL Java_org_phash_ni_audio_AudioHash_readAudio
(JNIEnv *env, jobject thisObj, jstring file, jint sr, jint nChannels, jfloat nbSeconds){

	const char *filename = env->GetStringUTFChars(file, NULL);

	int buflen;
	float *buf = ph_readaudio(filename, (int)sr, (int)nChannels, NULL, buflen, (float)nbSeconds);
	if (buf == NULL){
		env->ReleaseStringUTFChars(file, filename);
		return NULL;
	}
	jfloatArray arr = env->NewFloatArray((jsize)buflen);
	env->SetFloatArrayRegion(arr, 0, (jsize)buflen, (jfloat*)buf);
	env->ReleaseStringUTFChars(file, filename);
	return arr;
}

JNIEXPORT jintArray JNICALL Java_org_phash_ni_audio_AudioHash_audioHash
(JNIEnv *env, jobject thisObj, jfloatArray buf, jint sr){

	jfloat *ptrbuf = env->GetFloatArrayElements(buf, NULL);
	jsize  n = env->GetArrayLength(buf);

	int n_frames;
	uint32_t *ptrhash = ph_audiohash((float*)ptrbuf, (int)n, (int)sr, n_frames);
	if (ptrhash == NULL)
		return NULL;

	jintArray framesArr = env->NewIntArray((jsize)n_frames);
	env->SetIntArrayRegion(framesArr, 0, (jsize)n_frames, (jint*)ptrhash);

	free(framesArr);
	env->ReleaseFloatArrayElements(buf, ptrbuf, 0);
	
	return framesArr;
}


JNIEXPORT jdoubleArray JNICALL Java_org_phash_ni_audio_AudioHash_audioDistanceBER
(JNIEnv *env, jobject thisObj, jintArray hash1, jintArray hash2, jfloat threshold, jint blocksize){

	jint *phash1 = env->GetIntArrayElements(hash1, NULL);
	jsize n1 = env->GetArrayLength(hash1);
	jint *phash2 = env->GetIntArrayElements(hash2, NULL);
	jsize n2 = env->GetArrayLength(hash2);

	int N_cs;
	double *cs = ph_audio_distance_ber((uint32_t*)phash1, (int)n1, (uint32_t*)phash2, (int)n2,
									   (float)threshold, (int)blocksize, N_cs);

	
	jdoubleArray cscores = env->NewDoubleArray((jsize)N_cs);
	env->SetDoubleArrayRegion(cscores, 0, N_cs, cs);
	free(cs);
	
	return cscores;
}

#endif /* USE_AUDIO_HASH */


#ifdef USE_TEXT_HASH
	
JNIEXPORT jobjectArray JNICALL Java_org_phash_ni_text_TextHash_textHash
(JNIEnv *env, jobject thisObj, jstring file){
	const char *filename = env->GetStringUTFChars(file, NULL);

	int n;
	TxtHashPoint *txthash = ph_texthash(filename, &n);
	if (txthash == NULL){
		env->ReleaseStringUTFChars(file, filename);
		return NULL;
	}

	jobject txtHashPointObj = env->NewObject(txtHashPointClass, txtHashPointCtor, 0, 0); 
	jobjectArray arr = env->NewObjectArray((jsize)n, txtHashPointClass, txtHashPointObj);

	for (int i=0;i<n;i++){
		env->CallObjectMethod(txtHashPointObj, txtHashPoint_setDataMembers_longlong,
							  (jlong)txthash[i].hash, (jlong)txthash[i].index);
		env->SetObjectArrayElement(arr, (jsize)i, txtHashPointObj);
	}
	env->ReleaseStringUTFChars(file, filename);
	free(txthash);
	return arr;
}

JNIEXPORT jobjectArray JNICALL Java_org_phash_ni_text_TextHash_compareTextHashes
(JNIEnv *env, jobject thisObj, jobjectArray hash1, jobjectArray hash2){

	/* extract txt hashes from hash1 into c array */
	jsize n1 = env->GetArrayLength(hash1);
	TxtHashPoint *txthash1 = new TxtHashPoint[n1];
	for (int i=0;i<n1;i++){
		jobject point = env->GetObjectArrayElement(hash1, (jsize)i);
		txthash1[i].hash  = (long)env->CallObjectMethod(point, txtHashPoint_getHash_long);
		txthash1[i].index = (off_t)env->CallObjectMethod(point, txtHashPoint_getIndex_long);
	}
	
	jsize n2 = env->GetArrayLength(hash2);
	TxtHashPoint *txthash2 = new TxtHashPoint[n2];
	for (int i=0;i<n2;i++){
		jobject point = env->GetObjectArrayElement(hash2, (jsize)i);
		txthash2[i].hash = (long)env->CallObjectMethod(point, txtHashPoint_getHash_long);
		txthash2[i].index = (off_t)env->CallObjectMethod(point, txtHashPoint_getIndex_long);
	}
	
	/* compare */
	int n_matches;
	TxtMatch *matches = ph_compare_text_hashes(txthash1, n1, txthash2, n2, &n_matches);

	/* prepare jobjectarray to return to java */
	jobject txtMatchObj = env->NewObject(txtMatchClass, txtMatchCtor, 0, 0, 0);
	jobjectArray arr = env->NewObjectArray((jsize)n_matches, txtMatchClass, txtMatchObj);

	for (int i=0;i<n_matches;i++){
		env->CallObjectMethod(txtMatchObj, txtMatch_setDataMembers_longlonglong,
							  matches[i].first_index, matches[i].second_index, matches[i].length);
		env->SetObjectArrayElement(arr, (jsize)i, txtMatchObj);
		
	}

	delete [] txthash1;
	delete [] txthash2;
	free(matches);
	return arr;
}

#endif /* USE_TEXT_HASH */
	
#ifdef USE_VIDEO_HASH
	

JNIEXPORT jlongArray JNICALL Java_org_phash_ni_video_VideoHash_dctVideoHash
(JNIEnv *env, jobject thisObj, jstring file){
	const char *filename = env->GetStringUTFChars(file, NULL);

	int n;
	ulong64 *videohash = ph_dct_videohash(filename, n);

	jlongArray arr = env->NewLongArray((jsize)n);
	env->SetLongArrayRegion(arr, 0, (jsize)n, (jlong*)videohash);
	
	env->ReleaseStringUTFChars(file, filename);
	free(videohash);
	return arr;
}

JNIEXPORT jdouble JNICALL Java_org_phash_ni_video_VideoHash_dctVideoHashDistance
(JNIEnv *env, jobject thisObj, jlongArray hash1, jlongArray hash2, jint threshold){

	jlong *ptrhash1 = env->GetLongArrayElements(hash1, NULL);
	jsize n1 = env->GetArrayLength(hash1);
	
	jlong *ptrhash2 = env->GetLongArrayElements(hash2, NULL);
	jsize n2 = env->GetArrayLength(hash2);

	double d = ph_dct_videohash_dist((ulong64*)ptrhash1, n1, (ulong64*)ptrhash2, n2, (int)threshold);

	env->ReleaseLongArrayElements(hash1, ptrhash1, 0);
	env->ReleaseLongArrayElements(hash2, ptrhash2, 0);

	return (jdouble)d;
}

#endif /* USE_VIDEO_HASH */

