#include "config.h"
#include <jni.h>
#include "pHash-jni.h"
#include "pHash.h"
#include "pHash_MVPTree.h"

#ifdef HAVE_AUDIO_HASH
#include "audiophash.h"
#endif
 
jfieldID imHash_hash = NULL; 
jfieldID vidHash_hash = NULL; 
jfieldID audioHash_hash = NULL; 

jfieldID hash_filename = NULL; 

jclass imClass = NULL;
jclass vidClass = NULL;
jclass audioClass = NULL;

jmethodID imCtor = NULL;
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
	jclass cl;
	HashType hashType;
	hash_compareCB callback;
	jniHashType kind;
	jmethodID ctor;
	jfieldID hashID;
} jniHashes;


float image_distance(DP *pntA, DP *pntB)
{
    uint8_t htypeA = pntA->hash_type;
    uint8_t htypeB = pntB->hash_type;
    if (htypeA != htypeB)
        return -1.0;
    if (htypeA != UINT64ARRAY)
        return -1.0;
    if ((pntA->hash_length > 1) || (pntB->hash_length > 1))
        return -1.0;
    ulong64 *hashA = (ulong64*)pntA->hash;
    ulong64 *hashB = (ulong64*)pntB->hash;
    int res = ph_hamming_distance(*hashA, *hashB);
    return (float)res;
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
    for (int i=0;i<Nc;i++){
        if (ptrC[i] > maxC)
            maxC = ptrC[i];
    }
    double res = 1000*(1-maxC);
    return (float)res;
}


static jniHashes hashes[] = 
			{ 	
				{imClass, UINT64ARRAY, image_distance, IMAGE_HASH, imCtor, imHash_hash}, 
				{vidClass, UINT64ARRAY, image_distance, VIDEO_HASH, vidCtor, vidHash_hash}, 
				{audioClass, UINT32ARRAY, audio_distance, AUDIO_HASH, audioCtor, audioHash_hash},
			};

JNIEXPORT jboolean JNICALL Java_pHash_00024MVPTree_create
  (JNIEnv *e, jobject ob, jstring filename, jobjectArray hashArray)
{
	jint hashLen;
	if(filename == NULL || hashArray == NULL || (hashLen = e->GetArrayLength(hashArray)) == 0)
		return JNI_FALSE;

	const char *file = e->GetStringUTFChars(filename, 0);
	if(!file)
		return JNI_FALSE;

	MVPFile mvpfile;
	ph_mvp_init(&mvpfile);
	mvpfile.filename = file;
	jniHashType type;
	for(int i = 0; i < sizeof(hashes)/sizeof(hashes[0]); i++)
	{
		if(e->IsInstanceOf(hashArray, hashes[i].cl))
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
			e->ReleaseStringUTFChars(filename, file);
			return JNI_FALSE;
		}			
		jstring fname = (jstring)e->GetObjectField(hashObj, hash_filename);

		const char *path = e->GetStringUTFChars(fname, 0);
			
		hashlist[i]->id = strdup(path);
	
		switch(type)
		{
			case IMAGE_HASH:
				ulong64 tmphash;
				ph_dct_imagehash(path, tmphash);
				hashlist[i]->hash = (ulong64 *)malloc(sizeof(ulong64));
				if(!hashlist[i]->hash)
				{
					for(int i = 0; i < hashLen; i++)
					{
						if(hashlist[i])
							ph_free_datapoint(hashlist[i]);
					}

					free(hashlist);
					e->ReleaseStringUTFChars(filename, file);
					return JNI_FALSE;
				}
				*(ulong64 *)hashlist[i]->hash = tmphash;
				hashlist[i]->hash_length = 1;
				break;
			case VIDEO_HASH:
				ulong64 videoHash;
				ph_dct_videohash(path, videoHash);
                               	hashlist[i]->hash = (ulong64 *)malloc(sizeof(ulong64));
				if(!hashlist[i]->hash)
				{
					for(int i = 0; i < hashLen; i++)
					{
						if(hashlist[i])
							ph_free_datapoint(hashlist[i]);
					}

					free(hashlist);
					e->ReleaseStringUTFChars(filename, file);
					return JNI_FALSE;
				}
                                *(ulong64 *)hashlist[i]->hash = videoHash;
                                hashlist[i]->hash_length = 1;
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
					for(int i = 0; i < hashLen; i++)
					{
						if(hashlist[i])
							ph_free_datapoint(hashlist[i]);
					}

					free(hashlist);
					e->ReleaseStringUTFChars(filename, file);
					return JNI_FALSE;
				}
				break;
		}
		
		e->ReleaseStringUTFChars(fname, path);
	}

	for(int i = 0; i < hashLen; i++)
	{
		if(hashlist[i])
			ph_free_datapoint(hashlist[i]);
	}

	free(hashlist);
	e->ReleaseStringUTFChars(filename, file);
	return JNI_TRUE;

}

JNIEXPORT jobjectArray JNICALL Java_pHash_00024MVPTree_query
  (JNIEnv *e, jobject ob, jobject hashObj, jfloat radius, jint max)
{

	MVPFile mvpfile;
	jniHashType type;	
	ph_mvp_init(&mvpfile);
	
	jstring mvp = (jstring)e->GetObjectField(ob, e->GetFieldID(e->FindClass("pHash/MVPTree"), "mvpFile",
										"Ljava/lang/String;"));
	
	mvpfile.filename = e->GetStringUTFChars(mvp, 0);
	int i;
	for(i = 0; i < sizeof(hashes)/sizeof(hashes[0]); i++)
	{
		if(e->IsInstanceOf(hashObj, hashes[i].cl))
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
			query->hash_length = 1;
			ulong64 hash = (ulong64)e->GetLongField(hashObj, imHash_hash);
			query->hash = &hash;
			break;
		}
		case VIDEO_HASH:
		{
			query->hash_length = 1;
			ulong64 hash = (ulong64)e->GetLongField(hashObj, vidHash_hash);
                        query->hash = &hash;
			break;
		}
		case AUDIO_HASH:
		{
			hashList = (jintArray)e->GetObjectField(hashObj, audioHash_hash);
			query->hash_length = e->GetArrayLength(hashList);
			hash_list = e->GetIntArrayElements(hashList, NULL);
			query->hash = hash_list;
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
		jobject obj = e->NewObject(hashes[i].cl, hashes[i].ctor);
		ret = e->NewObjectArray(count, hashes[i].cl, obj);
		for(int j = 0; j < count; j++)
		{
			obj = e->GetObjectArrayElement(ret, j);
			jstring id = e->NewStringUTF(results[j]->id);
			e->SetObjectField(obj, hash_filename, id);
			switch(type)
			{
				case IMAGE_HASH:
					e->SetLongField(obj, hashes[i].hashID, *(jlong *)results[j]->hash);
					break;
				case VIDEO_HASH:
					e->SetLongField(obj, hashes[i].hashID, *(jlong *)results[j]->hash);
					break;
				case AUDIO_HASH:
					jintArray hashArray;
					e->SetIntArrayRegion(hashArray, 0, results[j]->hash_length, (jint *)results[j]->hash); 
					e->SetObjectField(obj, hashes[i].hashID, hashArray);
					break;
			}
		}
	}
	e->ReleaseStringUTFChars(mvp, mvpfile.filename);
	e->ReleaseStringUTFChars(hashStr, hash_file);
	ph_free_datapoint(query);
	for(int i = 0; i < max; i++)
	{
		if(results[i])
			ph_free_datapoint(results[i]);
	}
	free(results);
	return ret;
}

JNIEXPORT jboolean JNICALL Java_pHash_00024MVPTree_add
  (JNIEnv *e, jobject ob, jobjectArray hashArray)
{
	MVPFile mvpfile;
	jniHashType type;	
	ph_mvp_init(&mvpfile);

	if(hashArray == NULL)
		return JNI_FALSE;
	
	jstring mvp = (jstring)e->GetObjectField(ob, e->GetFieldID(e->FindClass("pHash/MVPTree"), "mvpFile",
										"Ljava/lang/String;"));
	
	mvpfile.filename = e->GetStringUTFChars(mvp, 0);
	int i;
	jobject hashObj = e->GetObjectArrayElement(hashArray, 0);
	for(i = 0; i < sizeof(hashes)/sizeof(hashes[0]); i++)
	{
		if(e->IsInstanceOf(hashObj, hashes[i].cl))
		{
			mvpfile.hashdist = hashes[i].callback;
			mvpfile.hash_type = hashes[i].hashType;
			type = hashes[i].kind;
			break;
		}
	}
	
	const char *hash_file = NULL;

	jsize len = e->GetArrayLength(hashArray);

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
			newHashes[j]->hash_length = 1;
			ulong64 hash = (ulong64)e->GetLongField(hashObj, imHash_hash);
			newHashes[j]->hash = (ulong64 *)malloc(sizeof(ulong64));
			*(ulong64 *)newHashes[j]->hash = hash;
			break;
		}
		case VIDEO_HASH:
		{
			newHashes[j]->hash_length = 1;
			ulong64 hash = (ulong64)e->GetLongField(hashObj, vidHash_hash);
                        newHashes[j]->hash = (ulong64 *)malloc(sizeof(ulong64));
			*(ulong64 *)newHashes[j]->hash = hash;
			break;
		}
		case AUDIO_HASH:
		{
			hashList = (jintArray)e->GetObjectField(hashObj, audioHash_hash);
			newHashes[j]->hash_length = e->GetArrayLength(hashList);
			hash_list = e->GetIntArrayElements(hashList, NULL);
			newHashes[j]->hash = hash_list;
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
	
	imClass = (jclass)e->NewGlobalRef(e->FindClass("ImageHash"));
	audioClass = (jclass)e->NewGlobalRef(e->FindClass("AudioHash"));
	vidClass = (jclass)e->NewGlobalRef(e->FindClass("VideoHash"));
        imHash_hash = e->GetFieldID(imClass, "hash", "J");
        audioHash_hash = e->GetFieldID(audioClass, "hash", "[I");
        vidHash_hash = e->GetFieldID(vidClass, "hash", "J");
	
	hash_filename = e->GetFieldID(e->FindClass("Hash"), "filename", "Ljava/lang/String;");

	imCtor = e->GetMethodID(imClass, "<init>", "()V");
	vidCtor = e->GetMethodID(vidClass, "<init>", "()V");
	audioCtor = e->GetMethodID(audioClass, "<init>", "()V");

}
JNIEXPORT void JNICALL Java_pHash_finalize
  (JNIEnv *e, jclass cl)
{
	e->DeleteGlobalRef(imClass);
	e->DeleteGlobalRef(vidClass);
	e->DeleteGlobalRef(audioClass);

}
JNIEXPORT jint JNICALL Java_pHash_imageDistance
  (JNIEnv *e, jclass cl, jobject hash1, jobject hash2)
{
	ulong64 imHash, imHash2;
	imHash = (ulong64)e->GetLongField(hash1, imHash_hash);
	imHash2 = (ulong64)e->GetLongField(hash2, imHash_hash);

	return ph_hamming_distance(imHash, imHash2);
	
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

JNIEXPORT jobject JNICALL Java_pHash_imageHash
  (JNIEnv *e, jclass cl, jstring f)
{
    
	const char *file = e->GetStringUTFChars(f,0);
    
    	ulong64 hash;
    	ph_dct_imagehash(file, hash);
	jobject imageHash = e->NewObject(imClass, imCtor);
	e->SetObjectField(imageHash, hash_filename, f);

	e->SetLongField(imageHash, imHash_hash, (jlong)hash);
    	e->ReleaseStringUTFChars(f,file);
	
	return imageHash;
	
}
#ifdef HAVE_VIDEO_HASH
JNIEXPORT jobject JNICALL Java_pHash_videoHash
  (JNIEnv *e, jclass cl, jstring f)
{
    
	const char *file = e->GetStringUTFChars(f,0);
    
	ulong64 hash;
	ph_dct_videohash(file, hash);

	jobject videoHash = e->NewObject(vidClass, vidCtor);
	e->SetObjectField(videoHash, hash_filename, f);

	e->SetLongField(videoHash, vidHash_hash, (jlong)hash);
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

	jintArray hashVals = (jintArray)e->GetObjectField(audioHash, audioHash_hash);
	
	e->SetIntArrayRegion(hashVals, 0, nbframes, (jint *)hash);
	free(hash);

	return audioHash;
	
}
#endif
