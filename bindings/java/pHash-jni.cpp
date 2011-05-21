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

#include "config.h"
#include <jni.h>
#include "pHash-jni.h"
#include "pHash.h"
#include "pHash_MVPTree.h"

#ifdef HAVE_AUDIO_HASH
#include "audiophash.h"
#endif
 
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

typedef enum ph_jni_hash_types
{
	IMAGE_HASH,
	VIDEO_HASH,
	AUDIO_HASH,
	TEXT_HASH
} jniHashType;

typedef struct ph_jni_hash_classes
{
	jclass *cl;
	HashType hashType;
	hash_compareCB callback;
	jniHashType kind;
	jmethodID *ctor;
	jfieldID *hashID;
} jniHashes;

float video_distance(DP *a, DP *b)
{
	ulong64 *hash1 = (ulong64 *)a->hash;
	ulong64 *hash2 = (ulong64 *)b->hash;
	double sim = ph_dct_videohash_dist(hash1, a->hash_length, hash2, b->hash_length, 21);
	return (float)sim;
}

float dctImage_distance(DP *pntA, DP *pntB)
{
	
	ulong64 *hashA = (ulong64*)pntA->hash;
	ulong64 *hashB = (ulong64*)pntB->hash;
	float res = ph_hamming_distance(*hashA, *hashB);
	return res;
}

float mhImage_distance(DP *a, DP *b)
{
	uint8_t *hashA = (uint8_t*)a->hash;
	uint8_t *hashB = (uint8_t*)b->hash;
	float res = ph_hammingdistance2(hashA,a->hash_length,hashB,b->hash_length)*10;
    	return exp(res);
}

float radialImage_distance(DP *a, DP *b)
{
	uint8_t *hashA = (uint8_t*)a->hash;
        uint8_t *hashB = (uint8_t*)b->hash;
	double pcc = -1;
	Digest x = {NULL,hashA, a->hash_length}, y = {NULL, hashB, b->hash_length};
        ph_crosscorr(x,y,pcc);
	return (float)pcc;

}
float audio_distance(DP *dpA, DP *dpB)
{
	uint32_t *hash1 = (uint32_t*)dpA->hash;
	int N1 = dpA->hash_length;
	uint32_t *hash2 = (uint32_t*)dpB->hash;
	int N2 = dpB->hash_length;

	float threshold = 0.30;
	int blocksize = 256;
	int Nc=0;
   	double *ptrC = ph_audio_distance_ber(hash1, N1, hash2, N2, threshold, blocksize, Nc);

    	double maxC = 0;
    	for (int i=0;i<Nc;i++)
	{
        	if (ptrC[i] > maxC)
            		maxC = ptrC[i];
    	}
	free(ptrC);
    	double res = 1000*(1-maxC);
    	return (float)res;
}


const static jniHashes hashes[] = 
			{ 	
				{&radialImClass, BYTEARRAY, radialImage_distance, IMAGE_HASH, &radialImCtor, &radialImHash_hash}, 
				{&mhImClass, BYTEARRAY, mhImage_distance, IMAGE_HASH, &mhImCtor, &mhImHash_hash}, 
				{&dctImClass, UINT64ARRAY, dctImage_distance, IMAGE_HASH, &dctImCtor, &dctImHash_hash},
				{&vidClass, UINT64ARRAY, video_distance, VIDEO_HASH, &vidCtor, &vidHash_hash}, 
				{&audioClass, UINT32ARRAY, audio_distance, AUDIO_HASH, &audioCtor, &audioHash_hash},
			};


JNIEXPORT jboolean JNICALL Java_org_phash_MVPTree_create
  (JNIEnv *e, jobject ob, jobjectArray hashArray)
{
	jint hashLen;
	if(hashArray == NULL || (hashLen = e->GetArrayLength(hashArray)) == 0)
		return JNI_FALSE;

	jstring mvp = (jstring)e->GetObjectField(ob, e->GetFieldID(e->FindClass("org/phash/MVPTree"), "mvpFile",
										"Ljava/lang/String;"));
	
	MVPFile mvpfile;
	ph_mvp_init(&mvpfile);
	mvpfile.filename = strdup(e->GetStringUTFChars(mvp, 0));
	jniHashType type;
	jobject htype = e->GetObjectArrayElement(hashArray, 0);
	for(int i = 0; i < sizeof(hashes)/sizeof(hashes[0]); i++)
	{
		if(e->IsInstanceOf(htype, *hashes[i].cl))
		{
			mvpfile.hashdist = hashes[i].callback;
			mvpfile.hash_type = hashes[i].hashType;
			type = hashes[i].kind;
			break;
		}
	}
	

	DP **hashlist = (DP **)malloc(hashLen*sizeof(DP*));
	for(jsize i = 0; i < hashLen; i++)
	{
		jobject hashObj = e->GetObjectArrayElement(hashArray, i);
		
		hashlist[i] = ph_malloc_datapoint(mvpfile.hash_type);
		if(!hashlist[i])
		{
			free(hashlist);
			e->ReleaseStringUTFChars(mvp, mvpfile.filename);
			return JNI_FALSE;
		}
		jstring fname = (jstring)e->GetObjectField(hashObj, hash_filename);

		const char *path = e->GetStringUTFChars(fname, 0);
			
		hashlist[i]->id = strdup(path);
	
		switch(type)
		{
			case IMAGE_HASH:
				if(e->IsInstanceOf(hashObj, dctImClass))
				{
					ulong64 tmphash;
					ph_dct_imagehash(path, tmphash);
					hashlist[i]->hash = (ulong64 *)malloc(sizeof(ulong64));
					*(ulong64 *)hashlist[i]->hash = tmphash;
					hashlist[i]->hash_length = 1;
				}
				else if(e->IsInstanceOf(hashObj, mhImClass))
				{
					int N;
					uint8_t *hash = ph_mh_imagehash(path, N);
					hashlist[i]->hash = hash;
					hashlist[i]->hash_length = N;
				}				
				else if(e->IsInstanceOf(hashObj, radialImClass))
				{
					Digest d;
					if(ph_image_digest(path, 1.0, 1.0, d, 180) >= 0)
					{
						hashlist[i]->hash = d.coeffs;
						hashlist[i]->hash_length = d.size;
					}
					 
				}
				
				break;
			case VIDEO_HASH:
				{
				int len;
				ulong64 *vhash = ph_dct_videohash(path, len);
				if(!vhash)
				{					
					for(int i = 0; i < hashLen; i++)
					{
						if(hashlist[i])
							ph_free_datapoint(hashlist[i]);
					}

					free(hashlist);
					e->ReleaseStringUTFChars(mvp, mvpfile.filename);
					return JNI_FALSE;
				}
                                hashlist[i]->hash = vhash;
                                hashlist[i]->hash_length = len;
				}
				break;
			case AUDIO_HASH:
				const float threshold = 0.30;
				const int block_size = 256;
				const int sr = 8000;
				const int channels = 1;
				int nbframes, N;
				float *buf = ph_readaudio(path,sr,channels,NULL,N);
				if(buf)
				{
					uint32_t *audioHash = ph_audiohash(buf,N,sr,nbframes);
					free(buf);
					hashlist[i]->hash_length = nbframes;
					hashlist[i]->hash = audioHash;
				} 
				else
				{
					free(hashlist);
					e->ReleaseStringUTFChars(mvp, mvpfile.filename);
					return JNI_FALSE;
				}
				break;
		}
		
		e->ReleaseStringUTFChars(fname, path);
	}



	MVPRetCode ret = ph_save_mvptree(&mvpfile, hashlist, hashLen);
	for(int i = 0; i < hashLen; i++)
	{
		if(hashlist[i])
		{
			free(hashlist[i]->hash);
			free(hashlist[i]->id);
		}
	}

	free(hashlist);
	e->ReleaseStringUTFChars(mvp, mvpfile.filename);
	free(mvpfile.filename);
	return (int)ret;
}

JNIEXPORT jobjectArray JNICALL Java_org_phash_MVPTree_query
  (JNIEnv *e, jobject ob, jobject hashObj, jfloat radius, 
	jfloat thresh, jint max)
{

	MVPFile mvpfile;
	jniHashType type;	
	ph_mvp_init(&mvpfile);
	
	jstring mvp = (jstring)e->GetObjectField(ob, e->GetFieldID(e->FindClass("org/phash/MVPTree"), "mvpFile",
										"Ljava/lang/String;"));
	
	mvpfile.filename = strdup(e->GetStringUTFChars(mvp, 0));
	int i;
	for(i = 0; i < sizeof(hashes)/sizeof(hashes[0]); i++)
	{
		if(e->IsInstanceOf(hashObj, *hashes[i].cl))
		{
			mvpfile.hashdist = hashes[i].callback;
			mvpfile.hash_type = hashes[i].hashType;
			type = hashes[i].kind;
			break;
		}
	}
	
	DP *query = ph_malloc_datapoint(mvpfile.hash_type);
	DP **results = (DP **)malloc(max*sizeof(DP *));
	const char *hash_file = NULL;
	jstring hashStr = (jstring)e->GetObjectField(hashObj, hash_filename);

	hash_file = e->GetStringUTFChars(hashStr, 0);

	query->id = strdup(hash_file);
	int count = 0;
	jint *hash_list = NULL;
	jintArray hashList = NULL;
	switch(type)
	{
		case IMAGE_HASH:
		{	
			if(e->IsInstanceOf(hashObj, dctImClass))
			{
				query->hash_length = 1;
				ulong64 hash = (ulong64)e->GetLongField(hashObj, dctImHash_hash);
				query->hash = &hash;
			}
			else if(e->IsInstanceOf(hashObj, mhImClass) )
			{
				jbyteArray hash = (jbyteArray)e->GetObjectField(hashObj, mhImHash_hash);
				query->hash_length = e->GetArrayLength(hash);
				jbyte *hashes = e->GetByteArrayElements(hash, NULL);
				query->hash = (uint8_t*)hashes;				
			}
			else if(e->IsInstanceOf(hashObj, radialImClass) )
			{
				jbyteArray hash = (jbyteArray)e->GetObjectField(hashObj, radialImHash_hash);
				query->hash_length = e->GetArrayLength(hash);
				jbyte *hashes = e->GetByteArrayElements(hash, NULL);
				query->hash = (uint8_t*)hashes;				
			}
			break;
		}
		case VIDEO_HASH:
		{
			jlongArray l = (jlongArray)e->GetObjectField(hashObj, vidHash_hash);
			query->hash_length = e->GetArrayLength(l);
			jlong *h = e->GetLongArrayElements(l, NULL);
                        query->hash = (ulong64*)h;
			break;
		}
		case AUDIO_HASH:
		{
			hashList = (jintArray)e->GetObjectField(hashObj, audioHash_hash);
			query->hash_length = e->GetArrayLength(hashList);
			hash_list = e->GetIntArrayElements(hashList, NULL);
			query->hash = (uint32_t*)hash_list;
			break;
		}
	}
	int res = ph_query_mvptree(&mvpfile, query, max, radius, 
					thresh, results, count);
	if(type == AUDIO_HASH)
		e->ReleaseIntArrayElements(hashList, hash_list, JNI_ABORT);
	jobjectArray ret;
	if(res != 0)
	{
		ret = NULL;
	}
	else
	{
		jobject iobj = e->NewObject(*hashes[i].cl, *hashes[i].ctor);
		ret = e->NewObjectArray(count, *hashes[i].cl, iobj);
		for(int j = 0; j < count; j++)
		{
			jobject obj = e->NewObject(*hashes[i].cl, *hashes[i].ctor);

			jstring id = e->NewStringUTF(results[j]->id);
			e->SetObjectField(obj, hash_filename, id);
			switch(type)
			{
				case IMAGE_HASH:
					if(e->IsInstanceOf(obj, dctImClass))
						e->SetLongField(obj, *hashes[i].hashID, *(jlong *)results[j]->hash);
					else if(e->IsInstanceOf(obj, mhImClass) ||
						e->IsInstanceOf(obj, radialImClass))
					{
						jbyteArray hash = e->NewByteArray(results[j]->hash_length);
						e->SetByteArrayRegion(hash, 0, results[j]->hash_length, (jbyte *)results[j]->hash);
						e->SetObjectField(obj, *hashes[i].hashID, hash);
					}	
				break;
				case VIDEO_HASH:
					e->SetLongField(obj, *hashes[i].hashID, *(jlong *)results[j]->hash);
					break;
				case AUDIO_HASH:
					jintArray hashArray = e->NewIntArray(results[j]->hash_length);
					e->SetIntArrayRegion(hashArray, 0, results[j]->hash_length, (jint *)results[j]->hash); 
					e->SetObjectField(obj, *hashes[i].hashID, hashArray);
					break;
			}
			e->SetObjectArrayElement(ret,j,obj);

		}

	}
	e->ReleaseStringUTFChars(mvp, mvpfile.filename);
	e->ReleaseStringUTFChars(hashStr, hash_file);
	ph_free_datapoint(query);
	for(int i = 0; i < count; i++)
	{
		if(results[i])
			ph_free_datapoint(results[i]);
	}
	free(results);
	free(mvpfile.filename);
	return ret;
}

JNIEXPORT jboolean JNICALL Java_org_phash_MVPTree_add
  (JNIEnv *e, jobject ob, jobjectArray hashArray)
{
	MVPFile mvpfile;
	jniHashType type;	
	ph_mvp_init(&mvpfile);
	jsize len;

	if(hashArray == NULL || (len = e->GetArrayLength(hashArray)) == 0)
		return JNI_FALSE;
	
	jstring mvp = (jstring)e->GetObjectField(ob, e->GetFieldID(e->FindClass("org/phash/MVPTree"), "mvpFile",
										"Ljava/lang/String;"));
	
	mvpfile.filename = strdup(e->GetStringUTFChars(mvp, 0));
	int i;
	jobject hashObj = e->GetObjectArrayElement(hashArray, 0);
	for(i = 0; i < sizeof(hashes)/sizeof(hashes[0]); i++)
	{
		if(e->IsInstanceOf(hashObj, *hashes[i].cl))
		{
			mvpfile.hashdist = hashes[i].callback;
			mvpfile.hash_type = hashes[i].hashType;
			type = hashes[i].kind;
			break;
		}
	}
	
	const char *hash_file = NULL;


	DP **newHashes = (DP **)malloc(len*sizeof(DP *));
	jintArray hashList = NULL;

	for(int j = 0; j < len; j++)
	{
		newHashes[j] = ph_malloc_datapoint(mvpfile.hash_type);
		hashObj = e->GetObjectArrayElement(hashArray, j);
		jstring hashStr = (jstring)e->GetObjectField(hashObj, hash_filename);

		hash_file = e->GetStringUTFChars(hashStr, 0);
		newHashes[j]->id = strdup(hash_file);
		e->ReleaseStringUTFChars(hashStr, hash_file);

	jint *hash_list = NULL;
	switch(type)
	{
		case IMAGE_HASH:
		{
			
			if(e->IsInstanceOf(hashObj, dctImClass))
			{
				newHashes[j]->hash_length = 1;
				ulong64 hash = (ulong64)e->GetLongField(hashObj, dctImHash_hash);
				newHashes[j]->hash = (ulong64 *)malloc(sizeof(ulong64));
				*(ulong64 *)newHashes[j]->hash = hash;
			}
			else if(e->IsInstanceOf(hashObj, mhImClass))
			{
				jbyteArray h = (jbyteArray)e->GetObjectField(hashObj, mhImHash_hash);
				newHashes[j]->hash_length = e->GetArrayLength(h);
				jbyte *hash = e->GetByteArrayElements(h, NULL);
				newHashes[j]->hash = (uint8_t*)hash;
			}
			else if(e->IsInstanceOf(hashObj, radialImClass))
			{
				jbyteArray h = (jbyteArray)e->GetObjectField(hashObj, radialImHash_hash);
				newHashes[j]->hash_length = e->GetArrayLength(h);
				jbyte *hash = e->GetByteArrayElements(h, NULL);
				newHashes[j]->hash = (uint8_t*)hash;
			}

			break;
		}
		case VIDEO_HASH:
		{
			jlongArray l = (jlongArray)e->GetObjectField(hashObj, vidHash_hash);
			newHashes[j]->hash_length = e->GetArrayLength(l);
			jlong *h = e->GetLongArrayElements(l, NULL);
                        newHashes[j]->hash = (ulong64*)h;
			break;
		}
		case AUDIO_HASH:
		{
			hashList = (jintArray)e->GetObjectField(hashObj, audioHash_hash);
			newHashes[j]->hash_length = e->GetArrayLength(hashList);
			hash_list = e->GetIntArrayElements(hashList, NULL);
			newHashes[j]->hash = (uint32_t*)hash_list;
			break;
		}
	}


	}

	int nbsaved = 0;
	int res = ph_add_mvptree(&mvpfile, newHashes, len, nbsaved);

	for(int j = 0; j < len; j++)
	{
		if(type == AUDIO_HASH)
		{
			hashObj = e->GetObjectArrayElement(hashArray, j);
			hashList = (jintArray)e->GetObjectField(hashObj, audioHash_hash);
			e->ReleaseIntArrayElements(hashList, (jint *)newHashes[j]->hash, JNI_ABORT);
		}
		ph_free_datapoint(newHashes[j]);
	}

	e->ReleaseStringUTFChars(mvp, mvpfile.filename);
	free(newHashes);
	free(mvpfile.filename);
	return JNI_TRUE;

}

JNIEXPORT void JNICALL Java_org_phash_pHash_pHashInit
  (JNIEnv *e, jclass cl)
{
	
	radialImClass = (jclass)e->NewGlobalRef(e->FindClass("org/phash/RadialImageHash"));
	dctImClass = (jclass)e->NewGlobalRef(e->FindClass("org/phash/DCTImageHash"));
	mhImClass = (jclass)e->NewGlobalRef(e->FindClass("org/phash/MHImageHash"));
	audioClass = (jclass)e->NewGlobalRef(e->FindClass("org/phash/AudioHash"));
	vidClass = (jclass)e->NewGlobalRef(e->FindClass("org/phash/VideoHash"));

        dctImHash_hash = e->GetFieldID(dctImClass, "hash", "J");
	radialImHash_hash = e->GetFieldID(radialImClass, "hash", "[B");
	mhImHash_hash = e->GetFieldID(mhImClass, "hash", "[B");
	audioHash_hash = e->GetFieldID(audioClass, "hash", "[I");
	vidHash_hash = e->GetFieldID(vidClass, "hash", "[J");

	hash_filename = e->GetFieldID(e->FindClass("org/phash/Hash"), "filename", "Ljava/lang/String;");

	dctImCtor = e->GetMethodID(dctImClass, "<init>", "()V");
	radialImCtor = e->GetMethodID(radialImClass, "<init>", "()V");
	mhImCtor = e->GetMethodID(mhImClass, "<init>", "()V");
	vidCtor = e->GetMethodID(vidClass, "<init>", "()V");
	audioCtor = e->GetMethodID(audioClass, "<init>", "()V");

}
JNIEXPORT void JNICALL Java_org_phash_pHash_cleanup
  (JNIEnv *e, jclass cl)
{
	e->DeleteGlobalRef(radialImClass);
	e->DeleteGlobalRef(mhImClass);
	e->DeleteGlobalRef(dctImClass);
	e->DeleteGlobalRef(vidClass);
	e->DeleteGlobalRef(audioClass);

}
JNIEXPORT jdouble JNICALL Java_org_phash_pHash_imageDistance
  (JNIEnv *e, jclass cl, jobject hash1, jobject hash2)
{
	if(e->IsInstanceOf(hash1, dctImClass) && e->IsInstanceOf(hash2, dctImClass))
	{
		ulong64 imHash, imHash2;
		imHash = (ulong64)e->GetLongField(hash1, dctImHash_hash);
		imHash2 = (ulong64)e->GetLongField(hash2, dctImHash_hash);

		return ph_hamming_distance(imHash, imHash2);
	}
	else if(e->IsInstanceOf(hash1, mhImClass) && e->IsInstanceOf(hash2, mhImClass))
	{
		jbyteArray h = (jbyteArray)e->GetObjectField(hash1, mhImHash_hash);
		jbyteArray h2 = (jbyteArray)e->GetObjectField(hash2, mhImHash_hash);
		int N = e->GetArrayLength(h);
		int N2 = e->GetArrayLength(h2);
		jbyte *hash = e->GetByteArrayElements(h, NULL);
		jbyte *hash2 = e->GetByteArrayElements(h2, NULL);
		double hd = ph_hammingdistance2((uint8_t*)hash, N, (uint8_t*)hash2, N2);
		e->ReleaseByteArrayElements(h, hash, 0);
		e->ReleaseByteArrayElements(h2, hash2, 0);
		return hd;	
	}
	else if(e->IsInstanceOf(hash1, radialImClass) && e->IsInstanceOf(hash2, radialImClass))
	{
		jbyteArray h = (jbyteArray)e->GetObjectField(hash1, radialImHash_hash);
		jbyteArray h2 = (jbyteArray)e->GetObjectField(hash2, radialImHash_hash);
		int N = e->GetArrayLength(h);
		int N2 = e->GetArrayLength(h2);
		jbyte *hash = e->GetByteArrayElements(h, NULL);
		jbyte *hash2 = e->GetByteArrayElements(h2, NULL);
		
		Digest x = {NULL,(uint8_t*)hash,N}, y = {NULL,(uint8_t*)hash2,N2};
		double pcc = -1;
		ph_crosscorr(x,y,pcc);
		e->ReleaseByteArrayElements(h, hash, 0);
		e->ReleaseByteArrayElements(h2, hash2, 0);
		return pcc;	
	}
	return -1;
}

#ifdef HAVE_AUDIO_HASH
JNIEXPORT jdouble JNICALL Java_org_phash_pHash_audioDistance
  (JNIEnv *e, jclass cl, jobject audioHash1, jobject audioHash2)
{

	if(audioHash1 == NULL || audioHash2 == NULL)
	{
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
#endif

JNIEXPORT jobject JNICALL Java_org_phash_pHash_dctImageHash
  (JNIEnv *e, jclass cl, jstring f)
{
    
	const char *file = e->GetStringUTFChars(f,0);
    
    	ulong64 hash;
    	ph_dct_imagehash(file, hash);
	jobject imageHash = e->NewObject(dctImClass, dctImCtor);
	e->SetObjectField(imageHash, hash_filename, f);

	e->SetLongField(imageHash, dctImHash_hash, (jlong)hash);
    	e->ReleaseStringUTFChars(f,file);
	
	return imageHash;
	
}

JNIEXPORT jobject JNICALL Java_org_phash_pHash_mhImageHash
  (JNIEnv *e, jclass cl, jstring f)
{
    
	const char *file = e->GetStringUTFChars(f,0);
    
	int N;
    	uint8_t *hash = ph_mh_imagehash(file, N);
	jobject imageHash = NULL;
	if(hash && N > 0)
	{
		imageHash = e->NewObject(mhImClass, mhImCtor);
		e->SetObjectField(imageHash, hash_filename, f);

		jbyteArray hashVals = e->NewByteArray(N);
		e->SetByteArrayRegion(hashVals, 0, N, (jbyte *)hash);
		e->SetObjectField(imageHash, mhImHash_hash, hashVals);
		free(hash);
	}
    	e->ReleaseStringUTFChars(f,file);
	
	return imageHash;
	
}
JNIEXPORT jobject JNICALL Java_org_phash_pHash_radialImageHash
  (JNIEnv *e, jclass cl, jstring f)
{
    
	const char *file = e->GetStringUTFChars(f,0);
    
	Digest dig;
	ph_image_digest(file, 1.0, 1.0, dig, 180);
	jobject imageHash = NULL;
	if(dig.coeffs && dig.size > 0)
	{
		imageHash = e->NewObject(radialImClass, radialImCtor);
		e->SetObjectField(imageHash, hash_filename, f);

		jbyteArray hashVals = e->NewByteArray(dig.size);
		e->SetByteArrayRegion(hashVals, 0, dig.size, (jbyte *)dig.coeffs);
		e->SetObjectField(imageHash, radialImHash_hash, hashVals);
		free(dig.coeffs);
	}
    	e->ReleaseStringUTFChars(f,file);
	
	return imageHash;
	
}
#ifdef HAVE_VIDEO_HASH

JNIEXPORT jdouble JNICALL Java_org_phash_pHash_videoDistance
  (JNIEnv *e, jclass cl, jobject vidHash1, jobject vidHash2, jint thresh)
{

	if(vidHash1 == NULL || vidHash2 == NULL)
	{
		return (jdouble)-1.0;
	}

	jint hash1_len, hash2_len;
	jlongArray hash1, hash2;
	hash1 = (jlongArray)e->GetObjectField(vidHash1, vidHash_hash);
	hash2 = (jlongArray)e->GetObjectField(vidHash2, vidHash_hash);
	hash1_len = e->GetArrayLength(hash1);
	hash2_len = e->GetArrayLength(hash2);
	if(hash1_len <= 0 || hash2_len <= 0)
	{
		return (jdouble)-1.0;
	}
	ulong64 *hash1_n, *hash2_n;
	
	hash1_n = (ulong64 *)e->GetLongArrayElements(hash1, 0);
	hash2_n = (ulong64 *)e->GetLongArrayElements(hash2, 0);
	jdouble sim = ph_dct_videohash_dist(hash1_n, hash1_len, hash2_n, hash2_len, thresh);
	e->ReleaseLongArrayElements(hash1, (jlong*)hash1_n, 0);
	e->ReleaseLongArrayElements(hash2, (jlong*)hash2_n, 0);
	return sim;	
}
JNIEXPORT jobject JNICALL Java_org_phash_pHash_videoHash
  (JNIEnv *e, jclass cl, jstring f)
{
    
	const char *file = e->GetStringUTFChars(f,0);
    	int len;
	ulong64 *hash =	ph_dct_videohash(file, len);

	jobject videoHash = e->NewObject(vidClass, vidCtor);
	e->SetObjectField(videoHash, hash_filename, f);

	jlongArray hashVals = e->NewLongArray(len);

	e->SetLongArrayRegion(hashVals, 0, len, (jlong *)hash);

	e->SetObjectField(videoHash, vidHash_hash, hashVals);
	free(hash);
    	e->ReleaseStringUTFChars(f,file);

	return videoHash;
	
}
#endif
#ifdef HAVE_AUDIO_HASH
JNIEXPORT jobject JNICALL Java_org_phash_pHash_audioHash
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
	buf = ph_readaudio(file,sr,channels,NULL,N); 
	if(!buf) 
	{
    		e->ReleaseStringUTFChars(f,file);
		return NULL;
	}
	hash = ph_audiohash(buf,N,sr,nbframes);
	if(!hash || nbframes <= 0) 
	{
    		free(buf);
		return NULL;
	}
	free(buf);

	jobject audioHash = e->NewObject(audioClass, audioCtor);
	e->SetObjectField(audioHash, hash_filename, f);
	e->ReleaseStringUTFChars(f,file);

	jintArray hashVals = e->NewIntArray(nbframes);

	e->SetIntArrayRegion(hashVals, 0, nbframes, (jint *)hash);
	e->SetObjectField(audioHash, audioHash_hash, hashVals);
	free(hash);

	return audioHash;
	
}
#endif
