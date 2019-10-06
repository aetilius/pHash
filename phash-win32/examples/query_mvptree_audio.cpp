#include <math.h>
#include <stdio.h>
#include "audiophash.h"
#include "pHash.h"

static int nb_calcs; /* metric to evaluate efficiency */
                     /* set to 0 before calling build, add, query */

float distfunc(DP *dpA, DP *dpB) {
    nb_calcs++;

    uint32_t *hash1 = (uint32_t *)dpA->hash;
    int N1 = dpA->hash_length;
    uint32_t *hash2 = (uint32_t *)dpB->hash;
    int N2 = dpB->hash_length;

    float threshold = 0.30f;
    int blocksize = 256;
    int Nc = 0;
    double *ptrC =
        ph_audio_distance_ber(hash1, N1, hash2, N2, threshold, blocksize, Nc);

    float maxC = 0;
    for (int i = 0; i < Nc; i++) {
        if (ptrC[i] > maxC) maxC = (float)ptrC[i];
    }
    if (ptrC) hfree(ptrC);

    double res = 1000 * (1 - maxC);
    return res;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        printf("not enough input args\n");
        printf(
            "usage: progname querydirectory dbname searchradius rnearest "
            "threshold nbsecs\n");
        return -1;
    }

    const char *dir_name =
        argv[1]; /* name of files in directory of query images */
    const char *filename = argv[2]; /* name of file to save db */

    const int sr = 8000;      /* sample to convert all files */
    const int nbchannels = 1; /* number channels to convert all files */

    MVPFile mvpfile; /* mvp configuration */
    mvpfile.filename = strdup(filename);
    mvpfile.hashdist = distfunc;
    mvpfile.hash_type = UINT32ARRAY;

    ULONG fl = 2;
    if (!HeapSetInformation(GetProcessHeap(), HeapEnableTerminationOnCorruption,
                            &fl, sizeof(ULONG))) {
        printf("unable to set heap information - %d\n", GetLastError());
        return -1;
    }

    /* get list of filenames in directory */
    int nbfiles = 0;
    printf("using db %s\n", filename);
    printf("using dir %s for query files\n", dir_name);
    char **files = ph_readfilenames(dir_name, nbfiles);
    if (!files) {
        printf("mem alloc error\n");
        return -1;
    }

    printf("nb query files = %d\n", nbfiles);

    /* search radius to traverse the tree */
    /* you will need to play around with this figure to find the optimal value
     */
    /* lower values will give you a more efficient search, but may not return
     * the result */
    /* high enough values should definitely return a result, if it is in there,
     * but might not be as efficient */
    float radius = 700.0; /* search radius to traverse the tree */
    if (argc >= 4) radius = atof(argv[3]);
    int rnearest = 1; /* max number results to retrieve - 1 for most efficient
                         searches (if each datapoint is unique */
    if (argc >= 5) rnearest = atoi(argv[4]);
    float threshold = 400.0; /*distance threshold for inclusion into list */
    if (argc >= 6) threshold = atof(argv[5]);
    float nbsecs = 45.0; /* number of seconds to read from audio file starting
                            from beginning  of track*/
    if (argc >= 7) nbsecs = atof(argv[6]);

    printf("radius %f, knearest %d, threshold %f, nb secs %f\n", radius,
           rnearest, threshold, nbsecs);

    DP *query = ph_malloc_datapoint(mvpfile.hash_type);

    DP **results = (DP **)malloc(rnearest * sizeof(DP *));
    if (results == NULL) {
        printf("unable to alloc memory\n");
        return -1;
    }

    float *sigbuf = (float *)malloc(1 << 21);
    int buflen = (1 << 21) / sizeof(float);

    uint32_t *hashbuf = (uint32_t *)malloc(1 << 14);
    int hashbuflength = (1 << 14) / sizeof(uint32_t);

    /* perform querys for each */
    int nbfound = 0, sum_calcs = 0, count = 0;
    for (int i = 0; i < nbfiles; i++) {
        printf("query file[%d]: %s\n", i, files[i]);
        int N = buflen;
        float *buf = ph_readaudio(files[i], sr, nbchannels, sigbuf, N, nbsecs);
        if (buf == NULL) {
            printf("unable to read audio\n");
            continue;
        }
        printf("buf length %d\n", N);

        int nbframes;
        uint32_t *hasht =
            ph_audiohash(buf, N, hashbuf, hashbuflength, sr, nbframes);
        if (hasht == NULL) {
            printf("unable to get hash\n");
            continue;
        }
        printf("hash %p to %p, length %d\n", hasht, hasht + nbframes, nbframes);

        query->id = files[i];
        query->hash = hasht;
        query->hash_length = nbframes;

        count++;
        nb_calcs = 0;
        nbfound = 0;
        printf("do query ...\n");
        MVPRetCode ret = ph_query_mvptree(&mvpfile, query, rnearest, radius,
                                          threshold, results, &nbfound);
        if (ret != PH_SUCCESS && ret != PH_ERRCAP) {
            printf("could not complete query - %d\n", ret);
            continue;
        }
        sum_calcs += nb_calcs;
        printf(" %d found, %d distance calcs\n", nbfound, nb_calcs);
        for (int j = 0; j < nbfound; j++) {
            printf("    %d  %s dist = %f\n", i, results[j]->id,
                   distfunc(results[j], query));
        }
        printf("******************************************\n");

        for (int j = 0; j < nbfound; j++) {
            hfree(results[j]->id);
            hfree(results[j]->hash);
            ph_free_datapoint(results[j]);
        }
    }
    float ave_calcs = (float)sum_calcs / (float)count;
    printf("ave calcs/query: %f\n", ave_calcs);

    ph_free_datapoint(query);
    free(results);

    hfree(sigbuf);
    hfree(hashbuf);

    for (int i = 0; i < nbfiles; i++) {
        sfree(files[i]);
    }
    hfree(files);
    sfree(mvpfile.filename);
    return 0;
}
