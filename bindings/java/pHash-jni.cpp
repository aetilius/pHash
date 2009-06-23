
#include "config.h"
#include "pHash-jni.h"
#include <jni.h>
#include "pHash.h"

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

typedef enum ph_jni_hash_types
{
	IMAGE_HASH,
	VIDEO_HASH,
	AUDIO_HASH,
	TEXT_HASH
} jniHashType;

typedef struct ph_jni_hash_classes
{
	jclass class;
	HashType hashType;
	hash_compareCB callback;
	jniHashType kind;
} jniHashes;

static jniHashes hashes[] = 
			{ 	
				{imClass, UINT64ARRAY, ph_hamming_distanceCB, IMAGE_HASH}, 
				{vidClass, UINT64ARRAY, ph_hamming_distanceCB, VIDEO_HASH}, 
				{audioClass, UINT32ARRAY, ph_audio_distance_berCB, AUDIO_HASH},
				{txtClass, BYTEARRAY, ph_text_distanceCB, TEXT_HASH}
			};

JNIEXPORT jboolean JNICALL Java_pHash_00024MVPTree_create
  (JNIEnv *e, jobject ob, jstring filename, jobjectArray hashArray)
{
	jint hashLen;
	if(filename == NULL || hashArray == NULL || (hashLen = e->GetArrayLength(hashArray)) == 0)
		return JNI_FALSE;

	const char *file = e->GetStringUTFChars(filename, 0);
	MVPFile mvpfile;
	ph_mvp_init(&mvpfile);
	mvpfile.filename = file;
	jniHashType type;
	for(int i = 0; i < sizeof(hashes)/sizeof(hashes[0]); i++)
	{
		if(e->IsInstanceOf(hashArray, hashes[i].class))
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
		jstring fname = e->GetObjectField(hashObj, hash_filename);

		const char *path = e->GetStringUTFChars(fname, 0);
		hashlist[i]->id = strdup(path);
	
		switch(type)
		{
			case IMAGE_HASH:
				ulong64 tmphash;
				ph_dct_imagehash(path, tmphash);
				hashlist[i]->hash = (ulong64 *)malloc(sizeof(ulong64));
				*(hashlist[i]->hash) = tmphash;
				hashlist[i]->hash_length = 1;
				break;
			case VIDEO_HASH:
				ulong64 videoHash;
				ph_dct_videohash(path, videoHash);
                               	hashlist[i]->hash = (ulong64 *)malloc(sizeof(ulong64));
                                *(hashlist[i]->hash) = videoHash;
                                hashlist[i]->hash_length = 1;
				break;
			case AUDIO_HASH:
				const float threshold = 0.30;
				const int block_size = 256;
				const int sr = 8000;
				const int channels = 1;
				int nbframes, N;
				float *buf = ph_read_audio(path,sr,channels,N);
				uint32_t *audioHash = ph_audiohash(buf,N,sr,nbframes);
				free(buf);
				haslist[i]->hash_length = nbframes;
				hashlist[i]->hash = audioHash;
				break;
			case TEXT_HASH:
				
				break;
		}
		
		e->ReleaseStringUTFChars(fname, path);

	}
	free(hashlist);
	e->ReleaseStringUTFChars(filename, file);	
}

JNIEXPORT jint JNICALL Java_pHash_pHashInit
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

JNIEXPORT jlong JNICALL Java_pHash_imageHash
  (JNIEnv *e, jclass cl, jstring f)
{
    
	const char *file = e->GetStringUTFChars(f,0);
    
    ulong64 hash;
    ph_dct_imagehash(file, hash);
    e->ReleaseStringUTFChars(f,file);

	return hash;
	
}
#ifdef HAVE_VIDEO_HASH
JNIEXPORT jlong JNICALL Java_pHash_videoHash
  (JNIEnv *e, jclass cl, jstring f)
{
    
	const char *file = e->GetStringUTFChars(f,0);
    
    ulong64 hash;
    ph_dct_videohash(file, hash);
    e->ReleaseStringUTFChars(f,file);

	return hash;
	
}
#endif
#ifdef HAVE_AUDIO_HASH
JNIEXPORT jintArray JNICALL Java_pHash_audioHash
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
    e->ReleaseStringUTFChars(f,file);
    hash = ph_audiohash(buf,N,sr,nbframes);
	if(!hash || nbframes <= 0) 
	{
    		free(buf);
		return NULL;
	}
	free(buf);

    ret = e->NewIntArray(nbframes);

    e->SetIntArrayRegion(ret, 0, nbframes, (jint *)hash);
    free(hash);

	return ret;
	
}
#endif
