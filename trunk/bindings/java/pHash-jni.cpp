#include "config.h"
#include <jni.h>
#include "pHash-jni.h"
#include "pHash.h"
#include "pHash_MVPTree.h"

#ifdef HAVE_AUDIO_HASH
#include "audiophash.h"
#endif
static jfieldID txtHash_hash = NULL; 
static jfieldID imHash_hash = NULL; 
static jfieldID vidHash_hash = NULL; 
static jfieldID audioHash_hash = NULL; 

static jfieldID txtHash_filename = NULL; 
static jfieldID imHash_filename = NULL; 
static jfieldID vidHash_filename = NULL; 
static jfieldID audioHash_filename = NULL; 
static jfieldID hash_filename = NULL; 

static jclass imClass = NULL;
static jclass vidClass = NULL;
static jclass audioClass = NULL;
static jclass txtClass = NULL;

static jmethodID imCtor = NULL;
static jmethodID vidCtor = NULL;
static jmethodID audioCtor = NULL;

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
				{imClass, UINT64ARRAY, image_distance, IMAGE_HASH}, 
				{vidClass, UINT64ARRAY, image_distance, VIDEO_HASH}, 
				{audioClass, UINT32ARRAY, audio_distance, AUDIO_HASH},
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

JNIEXPORT void JNICALL Java_pHash_pHashInit
  (JNIEnv *e, jclass cl)
{
	
	imClass = e->FindClass("pHash/ImageHash");
	txtClass = e->FindClass("pHash/TextHash");
	audioClass = e->FindClass("pHash/AudioHash");
	vidClass = e->FindClass("pHash/VideoHash");

        imHash_hash = e->GetFieldID(imClass, "hash", "J");
        audioHash_hash = e->GetFieldID(audioClass, "hash", "[I");
        txtHash_hash = e->GetFieldID(txtClass, "hash", "[I");
        vidHash_hash = e->GetFieldID(vidClass, "hash", "[J");
	
	hash_filename = e->GetFieldID(e->FindClass("pHash/Hash"), "filename", "Ljava/lang/String;");

	imCtor = e->GetMethodID(imClass, "<init>", "()V");
	vidCtor = e->GetMethodID(vidClass, "<init>", "()V");
	audioCtor = e->GetMethodID(audioClass, "<init>", "()V");

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
