#include <math.h>
#include <stdio.h>
#include "audiophash.h"
#include "pHash.h"

float distfunc(DP *dpA, DP *dpB) {
    uint32_t *hash1 = (uint32_t *)dpA->hash;
    int N1 = dpA->hash_length;
    uint32_t *hash2 = (uint32_t *)dpB->hash;
    int N2 = dpB->hash_length;

    float threshold = 0.30f;
    int blocksize = 256;
    int Nc = 0;
    double *ptrC =
        ph_audio_distance_ber(hash1, N1, hash2, N2, threshold, blocksize, Nc);

    double maxC = 0.0f;
    for (int i = 0; i < Nc; i++) {
        if (ptrC[i] > maxC) maxC = ptrC[i];
    }
    if (ptrC) free(ptrC);
    double res = 1000 * (1 - maxC);
    return (float)res;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        printf("not enough input args\n");
        printf("usage: %s directory dbname nbsecs\n", argv[0]);
        return -1;
    }
    const char *dir_name = argv[1]; /* name of dir to retrieve image files */
    const char *filename = argv[2]; /* name of file to save db */

    const int sr = 8000;      /* sample rate to convert all files */
    const int nbchannels = 1; /* number of channels to convert all files*/
    float nbsecs = 45.0f;     /* first number seconds of each file to read */
    if (argc >= 4) {
        nbsecs = atof(argv[3]);
    }

    MVPFile mvpfile; /* mvp tree configuration settings */
    mvpfile.filename = strdup(filename);
    mvpfile.hashdist = distfunc;
    mvpfile.hash_type = UINT32ARRAY;

    int nbfiles = 0;
    printf("dir name: %s\n", dir_name);
    char **files = ph_readfilenames(dir_name, nbfiles);
    if (!files) {
        printf("mem alloc error\n");
        return -1;
    }
    printf("number of files %d\n", nbfiles);
    DP **hashlist = (DP **)malloc(nbfiles * sizeof(DP *));
    if (hashlist == NULL) {
        printf("mem alloc error\n");
        return -1;
    }

    float *sigbuf = (float *)malloc(1 << 21);
    int buflen = (1 << 21) / sizeof(float);

    uint32_t *hashbuf = (uint32_t *)malloc(1 << 21);
    int hashbuflength = (1 << 21) / sizeof(uint32_t);
    uint32_t *hash = hashbuf;
    int hashbufleft = hashbuflength;

    int count = 0;
    for (int i = 0; i < nbfiles; i++) {
        printf("file[%d]: %s\n", i, files[i]);
        hashlist[count] = ph_malloc_datapoint(mvpfile.hash_type);
        int N = buflen;
        float *buf = ph_readaudio(files[i], sr, nbchannels, sigbuf, N, nbsecs);
        if (buf == NULL) {
            printf("cannot read file\n");
            continue;
        }
        int nbframes;
        uint32_t *hasht = ph_audiohash(buf, N, hash, hashbufleft, sr, nbframes);
        if (hasht == NULL) {
            printf("unable to get hash\n");
            continue;
        }
        hashlist[count]->id = files[i];
        hashlist[count]->hash = hash;
        hashlist[count]->hash_length = nbframes;
        count++;

        hash += nbframes;
        hashbufleft -= nbframes;
    }

    /* add each hash to list */
    printf("add %d files to file %s\n", count, filename);
    int nbsaved;
    MVPRetCode retcode = ph_add_mvptree(&mvpfile, hashlist, count, nbsaved);
    if (retcode != PH_SUCCESS) {
        printf("unable to add points to %s - %d\n", filename, retcode);
    }
    printf("save %d files out of %d\n", nbsaved, count);

    return 0;
}
