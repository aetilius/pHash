#ifndef _PHASH_JNI_H
#define _PHASH_JNI_H

#include <jni.h>

#cmakedefine USE_IMAGE_HASH
#cmakedefine USE_AUDIO_HASH
#cmakedefine USE_VIDEO_HASH
#cmakedefine USE_TEXT_HASH

#ifdef USE_IMAGE_HASH

/*
 * Class:     org_phash_ni_image_ImageHash
 * Method:    dctImageHash
 * Signature: (Ljava/lang/String;)J
 */
JNIEXPORT jlong JNICALL Java_org_phash_ni_image_ImageHash_dctImageHash
  (JNIEnv *, jobject, jstring);

/*
 * Class:     org_phash_ni_image_ImageHash
 * Method:    dctImageDistance
 * Signature: (JJ)I
 */
JNIEXPORT jint JNICALL Java_org_phash_ni_image_ImageHash_dctImageDistance
  (JNIEnv *, jobject, jlong, jlong);

/*
 * Class:     org_phash_ni_image_ImageHash
 * Method:    bmbImageHash
 * Signature: (Ljava/lang/String;)Lorg/phash/ni/image/BMBHash;
 */
JNIEXPORT jobject JNICALL Java_org_phash_ni_image_ImageHash_bmbImageHash
  (JNIEnv *, jobject, jstring);

/*
 * Class:     org_phash_ni_image_ImageHash
 * Method:    bmbDistance
 * Signature: (Lorg/phash/ni/image/BMBHash;Lorg/phash/ni/image/BMBHash;)I
 */
JNIEXPORT jint JNICALL Java_org_phash_ni_image_ImageHash_bmbDistance
  (JNIEnv *, jobject, jobject, jobject);

/*
 * Class:     org_phash_ni_image_ImageHash
 * Method:    radialFeatureVector
 * Signature: (Ljava/lang/String;IDD)Lorg/phash/ni/image/Features;
 */
JNIEXPORT jobject JNICALL Java_org_phash_ni_image_ImageHash_radialFeatureVector
  (JNIEnv *, jobject, jstring, jint, jdouble, jdouble);

/*
 * Class:     org_phash_ni_image_ImageHash
 * Method:    imageDigest
 * Signature: (Ljava/lang/String;DDI)Lorg/phash/ni/image/Digest;
 */
JNIEXPORT jobject JNICALL Java_org_phash_ni_image_ImageHash_imageDigest
  (JNIEnv *, jobject, jstring, jdouble, jdouble, jint);

/*
 * Class:     org_phash_ni_image_ImageHash
 * Method:    peakCrossCorr
 * Signature: (Lorg/phash/ni/image/Digest;Lorg/phash/ni/image/Digest;)D
 */
JNIEXPORT jdouble JNICALL Java_org_phash_ni_image_ImageHash_peakCrossCorr
  (JNIEnv *, jobject, jobject, jobject);

/*
 * Class:     org_phash_ni_image_ImageHash
 * Method:    mhImageHash
 * Signature: (Ljava/lang/String;FF)Lorg/phash/ni/image/MHash;
 */
JNIEXPORT jobject JNICALL Java_org_phash_ni_image_ImageHash_mhImageHash
  (JNIEnv *, jobject, jstring, jfloat, jfloat);

/*
 * Class:     org_phash_ni_image_ImageHash
 * Method:    mhImageDistance
 * Signature: (Lorg/phash/ni/image/MHash;Lorg/phash/ni/image/MHash;)I
 */
JNIEXPORT jint JNICALL Java_org_phash_ni_image_ImageHash_mhImageDistance
  (JNIEnv *, jobject, jobject, jobject);


#endif /* USE_IMAGE_HASH */

#ifdef USE_AUDIO_HASH
/*
 * Class:     org_phash_ni_audio_AudioHash
 * Method:    countAudioSamples
 * Signature: (Ljava/lang/String;II)I
 */
JNIEXPORT jint JNICALL Java_org_phash_ni_audio_AudioHash_countAudioSamples
  (JNIEnv *, jobject, jstring, jint, jint);

/*
 * Class:     org_phash_ni_audio_AudioHash
 * Method:    readAudio
 * Signature: (Ljava/lang/String;IIF)[F
 */
JNIEXPORT jfloatArray JNICALL Java_org_phash_ni_audio_AudioHash_readAudio
  (JNIEnv *, jobject, jstring, jint, jint, jfloat);

/*
 * Class:     org_phash_ni_audio_AudioHash
 * Method:    audioHash
 * Signature: ([FI)[I
 */
JNIEXPORT jintArray JNICALL Java_org_phash_ni_audio_AudioHash_audioHash
  (JNIEnv *, jobject, jfloatArray, jint);

/*
 * Class:     org_phash_ni_audio_AudioHash
 * Method:    audioDistanceBER
 * Signature: ([I[IFI)[D
 */
JNIEXPORT jdoubleArray JNICALL Java_org_phash_ni_audio_AudioHash_audioDistanceBER
  (JNIEnv *, jobject, jintArray, jintArray, jfloat, jint);

#endif /* USE_AUDIO_HASH */
	
#ifdef USE_TEXT_HASH
	
/*
 * Class:     org_phash_ni_text_TextHash
 * Method:    textHash
 * Signature: (Ljava/lang/String;)[Lorg/phash/ni/text/TxtHashPoint;
 */
JNIEXPORT jobjectArray JNICALL Java_org_phash_ni_text_TextHash_textHash
  (JNIEnv *, jobject, jstring);

/*
 * Class:     org_phash_ni_text_TextHash
 * Method:    compareTextHashes
 * Signature: ([Lorg/phash/ni/text/TxtHashPoint;[Lorg/phash/ni/text/TxtHashPoint;)[Lorg/phash/ni/text/TxtMatch;
 */
JNIEXPORT jobjectArray JNICALL Java_org_phash_ni_text_TextHash_compareTextHashes
  (JNIEnv *, jobject, jobjectArray, jobjectArray);

#endif /* USE_TEXT_HASH */

#ifdef USE_VIDEO_HASH
	
/*
 * Class:     org_phash_ni_video_VideoHash
 * Method:    dctVideoHash
 * Signature: (Ljava/lang/String;)[J
 */
JNIEXPORT jlongArray JNICALL Java_org_phash_ni_video_VideoHash_dctVideoHash
  (JNIEnv *, jobject, jstring);

/*
 * Class:     org_phash_ni_video_VideoHash
 * Method:    dctVideoHashDistance
 * Signature: ([J[JI)D
 */
JNIEXPORT jdouble JNICALL Java_org_phash_ni_video_VideoHash_dctVideoHashDistance
  (JNIEnv *, jobject, jlongArray, jlongArray, jint);

#endif /* USE_VIDEO_HASH */

#endif /* _PHASH_JNI_H */
