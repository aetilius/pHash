#include "config.h"
#include <jni.h>
#include "pHash-jni.h"
#include "pHash.h"
#include "pHash_MVPTree.h"
#include <assert.h>
#ifdef HAVE_AUDIO_HASH
#include "audiophash.h"
#endif
 
jfieldID dctImHash_hash = NULL; 
jfieldID mhImHash_hash = NULL; 
jfieldID vidHash_hash = NULL; 
jfieldID audioHash_hash = NULL; 

jfieldID hash_filename = NULL; 

jclass dctImClass = NULL;
jclass mhImClass = NULL;
jclass vidClass = NULL;
jclass audioClass = NULL;

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
	jfieldID hashID;
} jniHashes;

float video_distance(DP *a, DP *b)
{
	ulong64 *hash1 = (ulong64 *)a->hash;
	ulong64 *hash2 = (ulong64 *)b->hash;
	double sim = ph_dct_videohash_dist(hash1, a->hash_length, hash2, b->hash_length, 21);
	return (float)sim;
}

float image_distance(DP *pntA, DP *pntB)
{
	uint8_t htypeA = pntA->hash_type;
	uint8_t htypeB = pntB->hash_type;

	float res = 0;
	if (htypeA != htypeB)
        	return -1.0;
	if (htypeA != UINT64ARRAY && htypeA != BYTEARRAY)
        	return -1.0;
	if(htypeA == UINT64ARRAY)
	{
    		ulong64 *hashA = (ulong64*)pntA->hash;
    		ulong64 *hashB = (ulong64*)pntB->hash;
    		res = ph_hamming_distance(*hashA, *hashB);
	}
	else
	{
		uint8_t *hashA = (uint8_t*)pntA->hash;
		uint8_t *hashB = (uint8_t*)pntB->hash;
		res = ph_hammingdistance2(hashA,pntA->hash_length,hashB,pntB->hash_length)*1000;
	}
    	return res;
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


static jniHashes hashes[] = 
			{ 	
				{&mhImClass, BYTEARRAY, image_distance, IMAGE_HASH, &mhImCtor, mhImHash_hash}, 
				{&dctImClass, UINT64ARRAY, image_distance, IMAGE_HASH, &dctImCtor, dctImHash_hash}, 
				{&vidClass, UINT64ARRAY, video_distance, VIDEO_HASH, &vidCtor, vidHash_hash}, 
				{&audioClass, UINT32ARRAY, audio_distance, AUDIO_HASH, &audioCtor, audioHash_hash},
			};

JNIEXPORT jboolean JNICALL Java_MVPTree_create
  (JNIEnv *e, jobject ob, jobjectArray hashArray)
{
	jint hashLen;
	if(hashArray == NULL || (hashLen = e->GetArrayLength(hashArray)) == 0)
		return JNI_FALSE;

	jstring mvp = (jstring)e->GetObjectField(ob, e->GetFieldID(e->FindClass("MVPTree"), "mvpFile",
										"Ljava/lang/String;"));
	
	MVPFile mvpfile;
	ph_mvp_init(&mvpfile);
	mvpfile.filename = e->GetStringUTFChars(mvp, 0);
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
		
		hashlist[i] = ph_malloc_datapoint(mvpfile.hash_type, mvpfile.pathlength);
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
					assert(hash != NULL && N>0);			
					hashlist[i]->hash = hash;
					hashlist[i]->hash_length = N;
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
				float *buf = ph_readaudio(path,sr,channels,N);
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
			ph_free_datapoint(hashlist[i]);
	}

	free(hashlist);
	e->ReleaseStringUTFChars(mvp, mvpfile.filename);
	return ret == 0 ? JNI_TRUE : JNI_FALSE;
}

JNIEXPORT jobjectArray JNICALL Java_MVPTree_query
  (JNIEnv *e, jobject ob, jobject hashObj, jfloat radius, jint max)
{

	MVPFile mvpfile;
	jniHashType type;	
	ph_mvp_init(&mvpfile);
	
	jstring mvp = (jstring)e->GetObjectField(ob, e->GetFieldID(e->FindClass("MVPTree"), "mvpFile",
										"Ljava/lang/String;"));
	
	mvpfile.filename = e->GetStringUTFChars(mvp, 0);
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
	
	DP *query = ph_malloc_datapoint(mvpfile.hash_type, mvpfile.pathlength);
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
			else if(e->IsInstanceOf(hashObj, mhImClass))
			{
				jbyteArray hash = (jbyteArray)e->GetObjectField(hashObj, mhImHash_hash);
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
	int res = ph_query_mvptree(&mvpfile, query, max, radius, results, &count);
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

			printf("Result: %s\n", results[j]->id);
			jstring id = e->NewStringUTF(results[j]->id);
			e->SetObjectField(obj, hash_filename, id);
			switch(type)
			{
				case IMAGE_HASH:
					if(e->IsInstanceOf(obj, dctImClass))
						e->SetLongField(obj, hashes[i].hashID, *(jlong *)results[j]->hash);
					else if(e->IsInstanceOf(obj, mhImClass))
					{
						jbyteArray hash = e->NewByteArray(results[j]->hash_length);
						e->SetByteArrayRegion(hash, 0, results[j]->hash_length, (jbyte *)results[j]->hash);
						e->SetObjectField(obj, hashes[i].hashID, hash);
					}	
				break;
				case VIDEO_HASH:
					e->SetLongField(obj, hashes[i].hashID, *(jlong *)results[j]->hash);
					break;
				case AUDIO_HASH:
					jintArray hashArray = e->NewIntArray(results[j]->hash_length);
					e->SetIntArrayRegion(hashArray, 0, results[j]->hash_length, (jint *)results[j]->hash); 
					e->SetObjectField(obj, hashes[i].hashID, hashArray);
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
	return ret;
}

JNIEXPORT jboolean JNICALL Java_MVPTree_add
  (JNIEnv *e, jobject ob, jobjectArray hashArray)
{
	MVPFile mvpfile;
	jniHashType type;	
	ph_mvp_init(&mvpfile);
	jsize len;

	if(hashArray == NULL || (len = e->GetArrayLength(hashArray)) == 0)
		return JNI_FALSE;
	
	jstring mvp = (jstring)e->GetObjectField(ob, e->GetFieldID(e->FindClass("MVPTree"), "mvpFile",
										"Ljava/lang/String;"));
	
	mvpfile.filename = e->GetStringUTFChars(mvp, 0);
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
		newHashes[j] = ph_malloc_datapoint(mvpfile.hash_type, mvpfile.pathlength);
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

	int res = ph_add_mvptree(&mvpfile, newHashes, len);

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
	return JNI_TRUE;

}

JNIEXPORT void JNICALL Java_pHash_pHashInit
  (JNIEnv *e, jclass cl)
{
	
	dctImClass = (jclass)e->NewGlobalRef(e->FindClass("DCTImageHash"));
	mhImClass = (jclass)e->NewGlobalRef(e->FindClass("MHImageHash"));
	audioClass = (jclass)e->NewGlobalRef(e->FindClass("AudioHash"));
	vidClass = (jclass)e->NewGlobalRef(e->FindClass("VideoHash"));

        dctImHash_hash = e->GetFieldID(dctImClass, "hash", "J");
	mhImHash_hash = e->GetFieldID(mhImClass, "hash", "[B");
	audioHash_hash = e->GetFieldID(audioClass, "hash", "[I");
	vidHash_hash = e->GetFieldID(vidClass, "hash", "[J");

	hash_filename = e->GetFieldID(e->FindClass("Hash"), "filename", "Ljava/lang/String;");

	dctImCtor = e->GetMethodID(dctImClass, "<init>", "()V");
	mhImCtor = e->GetMethodID(mhImClass, "<init>", "()V");
	vidCtor = e->GetMethodID(vidClass, "<init>", "()V");
	audioCtor = e->GetMethodID(audioClass, "<init>", "()V");

}
JNIEXPORT void JNICALL Java_pHash_cleanup
  (JNIEnv *e, jclass cl)
{
	e->DeleteGlobalRef(mhImClass);
	e->DeleteGlobalRef(dctImClass);
	e->DeleteGlobalRef(vidClass);
	e->DeleteGlobalRef(audioClass);

}
JNIEXPORT jdouble JNICALL Java_pHash_imageDistance
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
	return -1;
}

#ifdef HAVE_AUDIO_HASH
JNIEXPORT jdouble JNICALL Java_pHash_audioDistance
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

JNIEXPORT jobject JNICALL Java_pHash_dctImageHash
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

JNIEXPORT jobject JNICALL Java_pHash_mhImageHash
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
#ifdef HAVE_VIDEO_HASH

JNIEXPORT jdouble JNICALL Java_pHash_videoDistance
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
JNIEXPORT jobject JNICALL Java_pHash_videoHash
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
JNIEXPORT jobject JNICALL Java_pHash_audioHash
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
