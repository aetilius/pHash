/*

 pHash, the open source perceptual hash library
 Copyright (C) 2010 Aetilius, Inc.
 All rights reserved.

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

 Evan Klinger - eklinger@phash.org
 David Starkweather - dstarkweather@phash.org

 */

#include <jni.h>
#include "pHash.h"
#include "pHashPro-jni.h"

jfieldID radialImHash_hash = NULL;
jfieldID dctImHash_hash = NULL;
jfieldID mhImHash_hash = NULL;
jfieldID vidHash_hash = NULL;
jfieldID audioHash_hash = NULL;

jfieldID hash_filename = NULL;

jclass radialImClass = NULL;
jclass dctImClass = NULL;
jclass mhImClass = NULL;
jclass vidClass = NULL;
jclass audioClass = NULL;

jmethodID radialImCtor = NULL;
jmethodID dctImCtor = NULL;
jmethodID mhImCtor = NULL;
jmethodID vidCtor = NULL;
jmethodID audioCtor = NULL;

typedef enum ph_jni_hash_types {
    IMAGE_HASH,
    VIDEO_HASH,
    AUDIO_HASH,
    TEXT_HASH
} jniHashType;

typedef struct ph_jni_hash_classes {
    jclass *cl;
    jniHashType kind;
    jmethodID *ctor;
    jfieldID *hashID;
} jniHashes;

const static jniHashes hashes[] = {
    {&radialImClass, IMAGE_HASH, &radialImCtor, &radialImHash_hash},
    {&mhImClass, IMAGE_HASH, &mhImCtor, &mhImHash_hash},
    {&dctImClass, IMAGE_HASH, &dctImCtor, &dctImHash_hash},
    {&vidClass, VIDEO_HASH, &vidCtor, &vidHash_hash},
    {&audioClass, AUDIO_HASH, &audioCtor, &audioHash_hash},
};

JNIEXPORT void JNICALL Java_org_phash_pHash_pHashInit(JNIEnv *e, jclass cl) {
    radialImClass =
        (jclass)e->NewGlobalRef(e->FindClass("org/phash/RadialImageHash"));
    dctImClass =
        (jclass)e->NewGlobalRef(e->FindClass("org/phash/DCTImageHash"));
    mhImClass = (jclass)e->NewGlobalRef(e->FindClass("org/phash/MHImageHash"));
    vidClass = (jclass)e->NewGlobalRef(e->FindClass("org/phash/VideoHash"));

    dctImHash_hash = e->GetFieldID(dctImClass, "hash", "J");
    radialImHash_hash = e->GetFieldID(radialImClass, "hash", "[B");
    mhImHash_hash = e->GetFieldID(mhImClass, "hash", "[B");
    audioHash_hash = e->GetFieldID(audioClass, "hash", "[I");
    vidHash_hash = e->GetFieldID(vidClass, "hash", "[J");

    hash_filename = e->GetFieldID(e->FindClass("org/phash/Hash"), "filename",
                                  "Ljava/lang/String;");

    dctImCtor = e->GetMethodID(dctImClass, "<init>", "()V");
    radialImCtor = e->GetMethodID(radialImClass, "<init>", "()V");
    mhImCtor = e->GetMethodID(mhImClass, "<init>", "()V");
    vidCtor = e->GetMethodID(vidClass, "<init>", "()V");
    audioCtor = e->GetMethodID(audioClass, "<init>", "()V");
}
JNIEXPORT void JNICALL Java_org_phash_pHash_cleanup(JNIEnv *e, jclass cl) {
    e->DeleteGlobalRef(radialImClass);
    e->DeleteGlobalRef(mhImClass);
    e->DeleteGlobalRef(dctImClass);
    e->DeleteGlobalRef(vidClass);
    e->DeleteGlobalRef(audioClass);
}
JNIEXPORT jdouble JNICALL Java_org_phash_pHash_imageDistance(JNIEnv *e,
                                                             jclass cl,
                                                             jobject hash1,
                                                             jobject hash2) {
    if (e->IsInstanceOf(hash1, dctImClass) &&
        e->IsInstanceOf(hash2, dctImClass)) {
        ulong64 imHash, imHash2;
        imHash = (ulong64)e->GetLongField(hash1, dctImHash_hash);
        imHash2 = (ulong64)e->GetLongField(hash2, dctImHash_hash);

        return ph_hamming_distance(imHash, imHash2);
    } else if (e->IsInstanceOf(hash1, mhImClass) &&
               e->IsInstanceOf(hash2, mhImClass)) {
        jbyteArray h = (jbyteArray)e->GetObjectField(hash1, mhImHash_hash);
        jbyteArray h2 = (jbyteArray)e->GetObjectField(hash2, mhImHash_hash);
        int N = e->GetArrayLength(h);
        int N2 = e->GetArrayLength(h2);
        jbyte *hash = e->GetByteArrayElements(h, NULL);
        jbyte *hash2 = e->GetByteArrayElements(h2, NULL);
        double hd =
            ph_hammingdistance2((uint8_t *)hash, N, (uint8_t *)hash2, N2);
        e->ReleaseByteArrayElements(h, hash, 0);
        e->ReleaseByteArrayElements(h2, hash2, 0);
        return hd;
    } else if (e->IsInstanceOf(hash1, radialImClass) &&
               e->IsInstanceOf(hash2, radialImClass)) {
        jbyteArray h = (jbyteArray)e->GetObjectField(hash1, radialImHash_hash);
        jbyteArray h2 = (jbyteArray)e->GetObjectField(hash2, radialImHash_hash);
        int N = e->GetArrayLength(h);
        int N2 = e->GetArrayLength(h2);
        jbyte *hash = e->GetByteArrayElements(h, NULL);
        jbyte *hash2 = e->GetByteArrayElements(h2, NULL);

        RadialHash x = {(uint8_t *)hash, N}, y = {(uint8_t *)hash2, N2};
        double pcc = -1;
        ph_peakcrosscorr(x, y, &pcc);
        e->ReleaseByteArrayElements(h, hash, 0);
        e->ReleaseByteArrayElements(h2, hash2, 0);
        return pcc;
    }
    return -1;
}

#ifdef HAVE_AUDIO_HASH
JNIEXPORT jdouble JNICALL Java_org_phash_pHash_audioDistance(
    JNIEnv *e, jclass cl, jobject audioHash1, jobject audioHash2) {
    if (audioHash1 == NULL || audioHash2 == NULL) {
        return (jdouble)-1.0;
    }
    const float threshold = 0.30;
    const int block_size = 256;
    int Nc;
    jdouble maxC = 0.0;
    double *pC;
    jint hash1_len, hash2_len;
    jintArray hash1, hash2;
    hash1 = (jintArray)e->GetObjectField(audioHash1, audioHash_hash);
    hash2 = (jintArray)e->GetObjectField(audioHash2, audioHash_hash);
    hash1_len = e->GetArrayLength(hash1);
    hash2_len = e->GetArrayLength(hash2);
    if (hash1_len <= 0 || hash2_len <= 0) {
        return (jdouble)-1.0;
    }
    uint32_t *hash1_n, *hash2_n;

    hash1_n = (uint32_t *)e->GetIntArrayElements(hash1, 0);
    hash2_n = (uint32_t *)e->GetIntArrayElements(hash2, 0);
    pC = ph_audio_distance_ber(hash1_n, hash1_len, hash2_n, hash2_len,
                               threshold, block_size, &Nc);
    e->ReleaseIntArrayElements(hash1, (jint *)hash1_n, 0);
    e->ReleaseIntArrayElements(hash2, (jint *)hash2_n, 0);
    maxC = 0.0;
    for (int j = 0; j < Nc; j++) {
        if (pC[j] > maxC) {
            maxC = pC[j];
        }
    }
    delete[] pC;
    return maxC;
}
#endif

JNIEXPORT jobject JNICALL Java_org_phash_pHash_dctImageHash(JNIEnv *e,
                                                            jclass cl,
                                                            jstring f) {
    const char *file = e->GetStringUTFChars(f, 0);

    ulong64 hash;
    int ret = ph_dct_imagehash(file, &hash);
    if (ret != 0) {
        e->ReleaseStringUTFChars(f, file);
        return NULL;
    }

    jobject imageHash = e->NewObject(dctImClass, dctImCtor);
    e->SetObjectField(imageHash, hash_filename, f);

    e->SetLongField(imageHash, dctImHash_hash, (jlong)hash);
    e->ReleaseStringUTFChars(f, file);

    return imageHash;
}

JNIEXPORT jobject JNICALL Java_org_phash_pHash_mhImageHash(JNIEnv *e, jclass cl,
                                                           jstring f) {
    const char *file = e->GetStringUTFChars(f, 0);

    int N;
    uint8_t *hash = ph_mh_imagehash(file, &N, 2, 1);
    jobject imageHash = NULL;
    if (hash && N > 0) {
        imageHash = e->NewObject(mhImClass, mhImCtor);
        e->SetObjectField(imageHash, hash_filename, f);

        jbyteArray hashVals = e->NewByteArray(N);
        e->SetByteArrayRegion(hashVals, 0, N, (jbyte *)hash);
        e->SetObjectField(imageHash, mhImHash_hash, hashVals);
        free(hash);
    }
    e->ReleaseStringUTFChars(f, file);

    return imageHash;
}
JNIEXPORT jobject JNICALL Java_org_phash_pHash_radialImageHash(JNIEnv *e,
                                                               jclass cl,
                                                               jstring f) {
    const char *file = e->GetStringUTFChars(f, 0);

    RadialHash rh;
    ph_radial_imagehash(file, &rh, 1.0, 180);
    jobject imageHash = NULL;
    if (rh.coeffs && rh.size > 0) {
        imageHash = e->NewObject(radialImClass, radialImCtor);
        e->SetObjectField(imageHash, hash_filename, f);

        jbyteArray hashVals = e->NewByteArray(rh.size);
        e->SetByteArrayRegion(hashVals, 0, rh.size, (jbyte *)rh.coeffs);
        e->SetObjectField(imageHash, radialImHash_hash, hashVals);
        free(rh.coeffs);
    }
    e->ReleaseStringUTFChars(f, file);

    return imageHash;
}
#ifdef HAVE_VIDEO_HASH

JNIEXPORT jdouble JNICALL Java_org_phash_pHash_videoDistance(
    JNIEnv *e, jclass cl, jobject vidHash1, jobject vidHash2, jint thresh) {
    if (vidHash1 == NULL || vidHash2 == NULL) {
        return (jdouble)-1.0;
    }

    jint hash1_len, hash2_len;
    jlongArray hash1, hash2;
    hash1 = (jlongArray)e->GetObjectField(vidHash1, vidHash_hash);
    hash2 = (jlongArray)e->GetObjectField(vidHash2, vidHash_hash);
    hash1_len = e->GetArrayLength(hash1);
    hash2_len = e->GetArrayLength(hash2);
    if (hash1_len <= 0 || hash2_len <= 0) {
        return (jdouble)-1.0;
    }
    ulong64 *hash1_n, *hash2_n;

    hash1_n = (ulong64 *)e->GetLongArrayElements(hash1, 0);
    hash2_n = (ulong64 *)e->GetLongArrayElements(hash2, 0);
    jdouble sim =
        ph_dct_videohash_dist(hash1_n, hash1_len, hash2_n, hash2_len, thresh);
    e->ReleaseLongArrayElements(hash1, (jlong *)hash1_n, 0);
    e->ReleaseLongArrayElements(hash2, (jlong *)hash2_n, 0);
    return sim;
}
JNIEXPORT jobject JNICALL Java_org_phash_pHash_videoHash(JNIEnv *e, jclass cl,
                                                         jstring f,
                                                         jdouble sigma) {
    const char *file = e->GetStringUTFChars(f, 0);
    int len;
    ulong64 *hashes = NULL;
    int ret = ph_dct_videohash_file(file, &hashes, &len, sigma);

    jobject videoHash = e->NewObject(vidClass, vidCtor);
    e->SetObjectField(videoHash, hash_filename, f);

    jlongArray hashVals = e->NewLongArray(len);

    e->SetLongArrayRegion(hashVals, 0, len, (jlong *)hashes);

    e->SetObjectField(videoHash, vidHash_hash, hashVals);
    free(hashes);
    e->ReleaseStringUTFChars(f, file);

    return videoHash;
}
#endif
#ifdef HAVE_AUDIO_HASH
JNIEXPORT jobject JNICALL Java_org_phash_pHash_audioHash(JNIEnv *e, jclass cl,
                                                         jstring f) {
    jintArray ret = NULL;
    const int sr = 8000;
    const int channels = 1;
    int N;
    int nbframes;
    float *buf = NULL;
    unsigned int *hash = NULL;
    const char *file = e->GetStringUTFChars(f, 0);
    buf = ph_readaudio(file, sr, channels, NULL, &N);
    if (!buf) {
        e->ReleaseStringUTFChars(f, file);
        return NULL;
    }
    hash = ph_audiohash(buf, N, sr, &nbframes);
    if (!hash || nbframes <= 0) {
        free(buf);
        return NULL;
    }
    free(buf);

    jobject audioHash = e->NewObject(audioClass, audioCtor);
    e->SetObjectField(audioHash, hash_filename, f);
    e->ReleaseStringUTFChars(f, file);

    jintArray hashVals = e->NewIntArray(nbframes);

    e->SetIntArrayRegion(hashVals, 0, nbframes, (jint *)hash);
    e->SetObjectField(audioHash, audioHash_hash, hashVals);
    free(hash);

    return audioHash;
}
#endif
