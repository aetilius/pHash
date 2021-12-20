/*

    pHash, the open source perceptual hash library
    Copyright (C) 2009 Aetilius, Inc.
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

#include "pHash.h"
#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <windows.h>
#include <cmath>
#include <utility>
#include "dirent.h"

#ifdef HAVE_VIDEO_HASH
#include "cimgffmpeg.h"
#endif

char phash_project[] = "%s. Copyright 2008-2009 Aetilius, Inc.";
char phash_version[255] = {0};

__declspec(dllexport) const char *ph_about() {
    if (phash_version[0] != 0) return phash_version;

    snprintf(phash_version, sizeof(phash_version), phash_project,
             PACKAGE_STRING);
    return phash_version;
}

#ifdef HAVE_IMAGE_HASH

__declspec(dllexport) int ph_radon_projections(const CImg<uint8_t> &img, int N,
                                               Projections &projs) {
    int width = img.width();
    int height = img.height();
    int D = (width > height) ? width : height;
    float x_center = (float)width / 2;
    float y_center = (float)height / 2;
    int x_off = std::round(x_center);
    int y_off = std::round(y_center);

    projs.R = new CImg<uint8_t>(N, D, 1, 1, 0);
    projs.nb_pix_perline = (int *)calloc(N, sizeof(int));

    if (!projs.R || !projs.nb_pix_perline) return EXIT_FAILURE;

    projs.size = N;

    CImg<uint8_t> *ptr_radon_map = projs.R;
    int *nb_per_line = projs.nb_pix_perline;

    for (int k = 0; k < N / 4 + 1; k++) {
        double theta = k * cimg::PI / N;
        double alpha = std::tan(theta);
        for (int x = 0; x < D; x++) {
            double y = alpha * (x - x_off);
            int yd = std::round(y);
            if ((yd + y_off >= 0) && (yd + y_off < height) && (x < width)) {
                *ptr_radon_map->data(k, x) = img(x, yd + y_off);
                nb_per_line[k] += 1;
            }
            if ((yd + x_off >= 0) && (yd + x_off < width) && (k != N / 4) &&
                (x < height)) {
                *ptr_radon_map->data(N / 2 - k, x) = img(yd + x_off, x);
                nb_per_line[N / 2 - k] += 1;
            }
        }
    }
    int j = 0;
    for (int k = 3 * N / 4; k < N; k++) {
        double theta = k * cimg::PI / N;
        double alpha = std::tan(theta);
        for (int x = 0; x < D; x++) {
            double y = alpha * (x - x_off);
            int yd = std::round(y);
            if ((yd + y_off >= 0) && (yd + y_off < height) && (x < width)) {
                *ptr_radon_map->data(k, x) = img(x, yd + y_off);
                nb_per_line[k] += 1;
            }
            if ((y_off - yd >= 0) && (y_off - yd < width) &&
                (2 * y_off - x >= 0) && (2 * y_off - x < height) &&
                (k != 3 * N / 4)) {
                *ptr_radon_map->data(k - j, x) =
                    img(-yd + y_off, -(x - y_off) + y_off);
                nb_per_line[k - j] += 1;
            }
        }
        j += 2;
    }

    return EXIT_SUCCESS;
}
__declspec(dllexport) int ph_feature_vector(const Projections &projs,
                                            Features &fv) {
    CImg<uint8_t> *ptr_map = projs.R;
    CImg<uint8_t> projection_map = *ptr_map;
    int *nb_perline = projs.nb_pix_perline;
    int N = projs.size;
    int D = projection_map.height();

    fv.features = (double *)malloc(N * sizeof(double));
    fv.size = N;
    if (!fv.features) return EXIT_FAILURE;

    double *feat_v = fv.features;
    double sum = 0.0;
    double sum_sqd = 0.0;
    for (int k = 0; k < N; k++) {
        double line_sum = 0.0;
        double line_sum_sqd = 0.0;
        int nb_pixels = nb_perline[k];
        for (int i = 0; i < D; i++) {
            line_sum += projection_map(k, i);
            line_sum_sqd += projection_map(k, i) * projection_map(k, i);
        }
        feat_v[k] = (line_sum_sqd / nb_pixels) -
                    (line_sum * line_sum) / (nb_pixels * nb_pixels);
        sum += feat_v[k];
        sum_sqd += feat_v[k] * feat_v[k];
    }
    double mean = sum / N;
    double var = sqrt((sum_sqd / N) - (sum * sum) / (N * N));

    for (int i = 0; i < N; i++) {
        feat_v[i] = (feat_v[i] - mean) / var;
    }

    return EXIT_SUCCESS;
}
__declspec(dllexport) int ph_dct(const Features &fv, Digest &digest) {
    int N = fv.size;
    const int nb_coeffs = 40;

    digest.coeffs = (uint8_t *)malloc(nb_coeffs * sizeof(uint8_t));
    if (!digest.coeffs) return EXIT_FAILURE;

    digest.size = nb_coeffs;

    double *R = fv.features;

    uint8_t *D = digest.coeffs;

    double D_temp[nb_coeffs];
    double max = 0.0;
    double min = 0.0;
    for (int k = 0; k < nb_coeffs; k++) {
        double sum = 0.0;
        for (int n = 0; n < N; n++) {
            double temp = R[n] * cos((cimg::PI * (2 * n + 1) * k) / (2 * N));
            sum += temp;
        }
        if (k == 0)
            D_temp[k] = sum / sqrt((double)N);
        else
            D_temp[k] = sum * SQRT_TWO / sqrt((double)N);
        if (D_temp[k] > max) max = D_temp[k];
        if (D_temp[k] < min) min = D_temp[k];
    }

    for (int i = 0; i < nb_coeffs; i++) {
        D[i] = (uint8_t)(UCHAR_MAX * (D_temp[i] - min) / (max - min));
    }

    return EXIT_SUCCESS;
}
__declspec(dllexport) int ph_crosscorr(const Digest &x, const Digest &y,
                                       double &pcc, double threshold) {
    int N = y.size;
    int result = 0;

    uint8_t *x_coeffs = x.coeffs;
    uint8_t *y_coeffs = y.coeffs;

    double *r = new double[N];
    double sumx = 0.0;
    double sumy = 0.0;
    for (int i = 0; i < N; i++) {
        sumx += x_coeffs[i];
        sumy += y_coeffs[i];
    }
    double meanx = sumx / N;
    double meany = sumy / N;
    double max = 0;
    for (int d = 0; d < N; d++) {
        double num = 0.0;
        double denx = 0.0;
        double deny = 0.0;
        for (int i = 0; i < N; i++) {
            num += (x_coeffs[i] - meanx) * (y_coeffs[(N + i - d) % N] - meany);
            denx += pow((x_coeffs[i] - meanx), 2);
            deny += pow((y_coeffs[(N + i - d) % N] - meany), 2);
        }
        r[d] = num / sqrt(denx * deny);
        if (r[d] > max) max = r[d];
    }
    delete[] r;
    pcc = max;
    if (max > threshold) result = 1;

    return result;
}

#ifdef max
#undef max
#endif

__declspec(dllexport) int ph_image_digest(const CImg<uint8_t> &img,
                                          double sigma, double gamma,
                                          Digest &digest, int N) {
    int result = EXIT_FAILURE;
    CImg<uint8_t> graysc;
    if (img.spectrum() >= 3) {
        graysc = img.get_RGBtoYCbCr().channel(0);
    } else if (img.spectrum() == 1) {
        graysc = img;
    } else {
        return result;
    }

    graysc.blur((float)sigma);

    (graysc / graysc.max()).pow(gamma);

    Projections projs;
    if (ph_radon_projections(graysc, N, projs) < 0) goto cleanup;

    Features features;
    if (ph_feature_vector(projs, features) < 0) goto cleanup;

    if (ph_dct(features, digest) < 0) goto cleanup;

    result = EXIT_SUCCESS;

cleanup:
    free(projs.nb_pix_perline);
    free(features.features);

    delete projs.R;
    return result;
}

#define max(a, b) (((a) > (b)) ? (a) : (b))

__declspec(dllexport) int ph_image_digest(const char *file, double sigma,
                                          double gamma, Digest &digest, int N) {
    CImg<uint8_t> *src = new CImg<uint8_t>(file);
    int result = ph_image_digest(*src, sigma, gamma, digest, N);
    delete src;
    return result;
}
__declspec(dllexport) int ph_compare_images(const CImg<uint8_t> &imA,
                                            const CImg<uint8_t> &imB,
                                            double &pcc, double sigma,
                                            double gamma, int N,
                                            double threshold) {
    int result = 0;
    Digest digestA;
    if (ph_image_digest(imA, sigma, gamma, digestA, N) < 0) goto cleanup;

    Digest digestB;
    if (ph_image_digest(imB, sigma, gamma, digestB, N) < 0) goto cleanup;

    if (ph_crosscorr(digestA, digestB, pcc, threshold) < 0) goto cleanup;

    if (pcc > threshold) result = 1;

cleanup:

    delete &imA;
    delete &imB;
    free(digestA.coeffs);
    free(digestB.coeffs);
    return result;
}
__declspec(dllexport) int ph_compare_images(const char *file1,
                                            const char *file2, double &pcc,
                                            double sigma, double gamma, int N,
                                            double threshold) {
    CImg<uint8_t> *imA = new CImg<uint8_t>(file1);
    CImg<uint8_t> *imB = new CImg<uint8_t>(file2);

    int res = ph_compare_images(*imA, *imB, pcc, sigma, gamma, N, threshold);

    return res;
}

__declspec(dllexport) CImg<float> *ph_dct_matrix(const int N) {
    CImg<float> *ptr_matrix = new CImg<float>(N, N, 1, 1, 1 / sqrt((float)N));
    float c1 = (float)sqrt(2.0 / N);
    for (int x = 0; x < N; x++) {
        for (int y = 1; y < N; y++) {
            *ptr_matrix->data(x, y) =
                c1 * (float)cos((cimg::PI / 2 / N) * y * (2 * x + 1));
        }
    }
    return ptr_matrix;
}

__declspec(dllexport) int ph_dct_imagehash(const char *file, ulong64 &hash) {
    if (!file) {
        return -1;
    }
    CImg<uint8_t> src;
    try {
        src.load(file);
    } catch (CImgIOException ex) {
        return -1;
    }
    CImg<float> meanfilter(7, 7, 1, 1, 1);
    CImg<float> img;
    if (src.spectrum() >= 3) {
        img = src.RGBtoYCbCr().channel(0).get_convolve(meanfilter);
    } else if (img.spectrum() == 1) {
        img = src.get_convolve(meanfilter);
    }

    img.resize(32, 32);
    CImg<float> *C = ph_dct_matrix(32);
    CImg<float> Ctransp = C->get_transpose();

    CImg<float> dctImage = (*C) * img * Ctransp;

    CImg<float> subsec = dctImage.crop(1, 1, 8, 8).unroll('x');
    ;

    float median = subsec.median();
    ulong64 one = 0x0000000000000001;
    hash = 0x0000000000000000;
    for (int i = 0; i < 64; i++) {
        float current = subsec(i);
        if (current > median) hash |= one;
        one = one << 1;
    }

    delete C;

    return 0;
}
#endif

#if defined(HAVE_VIDEO_HASH) && defined(HAVE_IMAGE_HASH)

__declspec(dllexport) ulong64 *ph_dct_videohash(const char *filename,
                                                int &Length) {
    CImgList<uint8_t> *keyframes = ph_getKeyFramesFromVideo(filename);
    if (keyframes == NULL) return NULL;

    Length = (int)keyframes->size();

    ulong64 *hash = (ulong64 *)malloc(sizeof(ulong64) * Length);
    CImg<float> *C = ph_dct_matrix(32);
    CImg<float> Ctransp = C->get_transpose();
    CImg<float> dctImage;
    CImg<float> subsec;
    CImg<uint8_t> currentframe;

    for (unsigned int i = 0; i < (unsigned int)keyframes->size(); i++) {
        currentframe = keyframes->at(i);
        currentframe.blur(1.0);
        dctImage = (*C) * (currentframe)*Ctransp;
        subsec = dctImage.crop(1, 1, 8, 8).unroll('x');
        float med = subsec.median();
        hash[i] = 0x0000000000000000;
        ulong64 one = 0x0000000000000001;
        for (int j = 0; j < 64; j++) {
            if (subsec[j] > med) hash[i] |= one;
            one = one << 1;
        }
    }

    keyframes->clear();
    delete keyframes;
    keyframes = NULL;
    delete C;
    C = NULL;
    return hash;
}

__declspec(dllexport) double ph_dct_videohash_dist(ulong64 *hashA, int N1,
                                                   ulong64 *hashB, int N2,
                                                   int threshold) {
    int den = (N1 <= N2) ? N1 : N2;
    int **C = new int *[N1 + 1];
    for (int i = 0; i < N1 + 1; i++) {
        C[i] = new int[N2 + 1];
    }

    for (int i = 0; i < N1 + 1; i++) {
        C[i][0] = 0;
    }
    for (int j = 0; j < N2 + 1; j++) {
        C[0][j] = 0;
    }
    for (int i = 1; i < N1 + 1; i++) {
        for (int j = 1; j < N2 + 1; j++) {
            int d = ph_hamming_distance(hashA[i - 1], hashB[j - 1]);
            if (d <= threshold) {
                C[i][j] = C[i - 1][j - 1] + 1;
            } else {
                C[i][j] =
                    ((C[i - 1][j] >= C[i][j - 1])) ? C[i - 1][j] : C[i][j - 1];
            }
        }
    }

    double result = (double)(C[N1][N2]) / (double)(den);

    for (int i = 0; i < N1 + 1; i++) {
        delete C[i];
    }
    delete[] C;

    return result;
}

#endif

__declspec(dllexport) int ph_hamming_distance(const ulong64 hash1,
                                              const ulong64 hash2) {
    ulong64 x = hash1 ^ hash2;
    const ulong64 m1 = 0x5555555555555555ULL;
    const ulong64 m2 = 0x3333333333333333ULL;
    const ulong64 h01 = 0x0101010101010101ULL;
    const ulong64 m4 = 0x0f0f0f0f0f0f0f0fULL;
    x -= (x >> 1) & m1;
    x = (x & m2) + ((x >> 2) & m2);
    x = (x + (x >> 4)) & m4;
    return int((x * h01) >> 56);
}

__declspec(dllexport) DP *ph_malloc_datapoint(HashType type) {
    DP *dp = (DP *)malloc(sizeof(DP));
    dp->hash = NULL;
    dp->id = NULL;
    dp->path = NULL;
    dp->hash_type = type;
    return dp;
}
__declspec(dllexport) void ph_free_datapoint(DP *dp) {
    if (!dp) return;
    hfree(dp);
    return;
}

__declspec(dllexport) DP **ph_read_imagehashes(const char *dirname,
                                               int pathlength, int &count) {
    count = 0;
    struct dirent *dir_entry;
    DIR *dir = opendir(dirname);
    if (!dir) exit(1);

    while ((dir_entry = readdir(dir)) != 0) {
        if (strcmp(dir_entry->d_name, ".") && strcmp(dir_entry->d_name, "..")) {
            count++;
        }
    }

    DP **hashlist = (DP **)malloc(count * sizeof(DP **));
    if (!hashlist) exit(1);

    DP *dp = NULL;
    int index = 0;
    errno = 0;
    ulong64 tmphash = 0;
    char path[100];
    path[0] = '\0';
    rewinddir(dir);
    while ((dir_entry = readdir(dir)) != 0) {
        if (strcmp(dir_entry->d_name, ".") && strcmp(dir_entry->d_name, "..")) {
            strcat(path, dirname);
            strcat(path, "/");
            strcat(path, dir_entry->d_name);
            if (ph_dct_imagehash(path, tmphash) < 0)  // calculate the hash
                continue;
            dp = ph_malloc_datapoint(UINT64ARRAY);
            dp->id = strdup(path);
            dp->hash = (void *)&tmphash;
            dp->hash_length = 1;
            hashlist[index++] = dp;
        }
        errno = 0;
        path[0] = '\0';
    }
    if (errno) exit(1);
    closedir(dir);
    return hashlist;
}

__declspec(dllexport) char **ph_readfilenames(const char *dirname, int &count) {
    count = 0;
    struct dirent *dir_entry;
    DIR *dir = opendir(dirname);
    if (!dir) return NULL;

    /*count files */
    while ((dir_entry = readdir(dir)) != NULL) {
        if (strcmp(dir_entry->d_name, ".") && strcmp(dir_entry->d_name, ".."))
            count++;
    }

    /* alloc list of files */
    char **files = (char **)malloc(count * sizeof(char *));
    if (!files) {
        return NULL;
    }

    int index = 0;
    char path[512];
    path[0] = '\0';
    rewinddir(dir);
    while ((dir_entry = readdir(dir)) != 0) {
        if (strcmp(dir_entry->d_name, ".") && strcmp(dir_entry->d_name, "..")) {
            snprintf(path, sizeof(path), "%s\\%s", dirname, dir_entry->d_name);
            files[index++] = strdup(path);
        }
    }
    count = index;
    closedir(dir);
    return files;
}

CImg<float> *GetMHKernel(float alpha, float level) {
    int sigma = (int)4 * pow((float)alpha, (float)level);
    static CImg<float> *pkernel = NULL;
    float xpos, ypos, A;
    if (!pkernel) {
        pkernel = new CImg<float>(2 * sigma + 1, 2 * sigma + 1, 1, 1, 0);
        cimg_forXY(*pkernel, X, Y) {
            xpos = pow((float)alpha, (float)-level) * (X - sigma);
            ypos = pow((float)alpha, (float)-level) * (Y - sigma);
            A = xpos * xpos + ypos * ypos;
            pkernel->atXY(X, Y) = (2 - A) * exp(-A / 2);
        }
    }
    return pkernel;
}
__declspec(dllexport) uint8_t *ph_mh_imagehash(const char *filename, int &N,
                                               float alpha, float lvl) {
    if (filename == NULL) {
        return NULL;
    }
    uint8_t *hash = (uint8_t *)malloc(72 * sizeof(uint8_t));
    if (!hash) {
        return NULL;
    }
    N = 72;

    CImg<uint8_t> src(filename);
    CImg<uint8_t> img;
    if (img.spectrum() == 3) {
        img = src.RGBtoYCbCr()
                  .channel(0)
                  .blur(1.0)
                  .resize(512, 512, 1, 1, 5)
                  .get_equalize(256);
    } else {
        img = src.channel(0)
                  .blur(1.0)
                  .resize(512, 512, 1, 1, 5)
                  .get_equalize(256);
    }
    src.clear();

    CImg<float> *pMHKernel = GetMHKernel(alpha, lvl);

    CImg<float> fresp = img.get_correlate(*pMHKernel);
    img.clear();
    fresp.normalize(0, 1.0);
    CImg<float> blocks(31, 31, 1, 1, 0);
    for (int rindex = 0; rindex < 31; rindex++) {
        for (int cindex = 0; cindex < 31; cindex++) {
            blocks(rindex, cindex) =
                fresp
                    .get_crop(rindex * 16, cindex * 16, rindex * 16 + 16 - 1,
                              cindex * 16 + 16 - 1)
                    .sum();
        }
    }

    int hash_index;
    int nb_ones = 0, nb_zeros = 0;
    int bit_index = 0;
    unsigned char hashbyte = 0;
    for (int rindex = 0; rindex < 31 - 2; rindex += 4) {
        CImg<float> subsec;
        for (int cindex = 0; cindex < 31 - 2; cindex += 4) {
            subsec = blocks.get_crop(cindex, rindex, cindex + 2, rindex + 2)
                         .unroll('x');
            float ave = subsec.mean();
            cimg_forX(subsec, I) {
                hashbyte <<= 1;
                if (subsec(I) > ave) {
                    hashbyte |= 0x01;
                    nb_ones++;
                } else {
                    nb_zeros++;
                }
                bit_index++;
                if ((bit_index % 8) == 0) {
                    hash_index = (int)(bit_index / 8) - 1;
                    hash[hash_index] = hashbyte;
                    hashbyte = 0x00;
                }
            }
        }
    }
    return hash;
}

__declspec(dllexport) int ph_bitcount8(uint8_t val) {
    int num = 0;
    /*
    while (val){
        num += val & 1;
        val >>= 1;
    }
    */
    while (val) {
        num++;
        val &= val - 1;
    }
    return num;
}

__declspec(dllexport) double ph_hammingdistance2(uint8_t *hashA, int lenA,
                                                 uint8_t *hashB, int lenB) {
    if (lenA != lenB) {
        return -1.0;
    }
    if ((hashA == NULL) || (hashB == NULL) || (lenA <= 0)) {
        return -1.0;
    }
    double dist = 0;
    uint8_t D = 0;
    for (int i = 0; i < lenA; i++) {
        D = hashA[i] ^ hashB[i];
        dist = dist + (double)ph_bitcount8(D);
    }
    double bits = (double)lenA * 8;
    return dist / bits;
}

int ph_selectvantagepoints(MVPFile *m, DP **points, int N, int &sv1_pos,
                           int &sv2_pos) {
    if (N > 2) {
        sv1_pos = 0;
        sv2_pos = 1;
        float maxdist = 0.0f;
        float d;
        for (int i = 0; i < N; i++) { /* find 2 points furthest apart */
            for (int j = i + 1; j < N; j++) {
                d = m->hashdist(points[i], points[j]);
                if (d > maxdist) {
                    maxdist = d;
                    sv1_pos = i;
                    sv2_pos = j;
                }
            }
        }
    } else if (N > 1) {
        sv1_pos = 0;
        sv2_pos = 1;
    } else if (N > 0) {
        sv1_pos = 0;
        sv2_pos = -1;
    } else {
        sv1_pos = -1;
        sv2_pos = -1;
    }
    return 0;
}

__declspec(dllexport) DP *ph_read_datapoint(MVPFile *m) {
    DP *dp = NULL;
    uint8_t active;
    uint16_t byte_len;
    uint16_t id_len;
    uint32_t hash_len;
    int type = m->hash_type;
    int PathLength = m->pathlength;
    off_t file_pos = m->file_pos;
    off_t offset_mask = m->pgsize - 1;
    memcpy(&active, &(m->buf[file_pos & offset_mask]), sizeof(uint8_t));
    file_pos++;
    memcpy(&byte_len, &(m->buf[file_pos & offset_mask]), sizeof(uint16_t));
    file_pos += sizeof(uint16_t);
    if ((active == 0) || (byte_len == 0)) {
        return dp;
    }

    dp = ph_malloc_datapoint(m->hash_type);
    dp->path = (float *)malloc(PathLength * sizeof(float));

    memcpy(&id_len, &(m->buf[file_pos & offset_mask]), sizeof(uint16_t));
    file_pos += sizeof(uint16_t);
    dp->id = (char *)malloc((id_len + 1) * sizeof(uint8_t));
    memcpy(dp->id, &(m->buf[file_pos & offset_mask]), id_len);
    dp->id[id_len] = '\0';
    file_pos += id_len;

    memcpy(&hash_len, &(m->buf[file_pos & offset_mask]), sizeof(uint32_t));
    dp->hash_length = hash_len;
    file_pos += sizeof(uint32_t);

    dp->hash = malloc(hash_len * type);
    memcpy(dp->hash, &(m->buf[file_pos & offset_mask]), hash_len * type);
    file_pos += hash_len * type;
    memcpy(dp->path, &(m->buf[file_pos & offset_mask]),
           PathLength * sizeof(float));
    file_pos += PathLength * sizeof(float);

    m->file_pos = file_pos;

    return dp;
}

__declspec(dllexport) int ph_sizeof_dp(DP *dp, MVPFile *m) {
    if (!dp) return (sizeof(uint8_t) + sizeof(uint16_t));
    return (int(strlen(dp->id)) + 9 + int((dp->hash_length) * (m->hash_type)) +
            int((m->pathlength) * sizeof(float)));
}

__declspec(dllexport) off_t ph_save_datapoint(DP *dp, MVPFile *m) {
    uint8_t active = 1;
    uint16_t byte_len = 0;
    off_t point_pos = m->file_pos;
    off_t offset_mask = m->pgsize - 1;
    if (dp == NULL) {
        active = 0;
        memcpy(&(m->buf[m->file_pos & offset_mask]), &active, 1);
        m->file_pos++;
        memcpy(&(m->buf[m->file_pos & offset_mask]), &byte_len,
               sizeof(uint16_t));
        m->file_pos += sizeof(uint16_t);
        return point_pos;
    }
    int type = m->hash_type;
    int PathLength = m->pathlength;
    uint16_t id_len = (uint16_t)strlen(dp->id);
    uint32_t hash_len = dp->hash_length;
    byte_len = (uint16_t)(id_len + 6 + hash_len * (m->hash_type) +
                          PathLength * sizeof(float));

    memcpy(&(m->buf[m->file_pos & offset_mask]), &active, 1);
    m->file_pos++;

    memcpy(&(m->buf[m->file_pos & offset_mask]), &byte_len, sizeof(uint16_t));
    m->file_pos += sizeof(uint16_t);

    memcpy(&(m->buf[m->file_pos & offset_mask]), &id_len, sizeof(uint16_t));
    m->file_pos += sizeof(uint16_t);

    memcpy(&(m->buf[m->file_pos & offset_mask]), (dp->id), id_len);
    m->file_pos += id_len;

    memcpy(&(m->buf[m->file_pos & offset_mask]), &hash_len, sizeof(uint32_t));
    m->file_pos += sizeof(uint32_t);

    memcpy(&(m->buf[m->file_pos & offset_mask]), (dp->hash), hash_len * type);
    m->file_pos += hash_len * type;

    memcpy(&(m->buf[m->file_pos & offset_mask]), (dp->path),
           sizeof(float) * PathLength);
    m->file_pos += PathLength * sizeof(float);

    return point_pos;
}

/* getpagesize for windows */
DWORD getpagesize(void) {
    static long g_pagesize = 0;
    if (!g_pagesize) {
        SYSTEM_INFO system_info;
        GetSystemInfo(&system_info);
        g_pagesize = system_info.dwPageSize;
    }
    return g_pagesize;
}
/* get alloc size for MapViewOfFile function */
DWORD getregionsize(void) {
    static long g_regionsize = 0;
    if (!g_regionsize) {
        SYSTEM_INFO system_info;
        GetSystemInfo(&system_info);
        g_regionsize = system_info.dwAllocationGranularity;
    }
    return g_regionsize;
}

MVPRetCode _ph_map_mvpfile(uint8_t filenumber, off_t offset, MVPFile *m,
                           MVPFile *m2, int use_existing) {
    char objname[256];
    char *fm_objname = NULL;
    m2->filename = strdup(m->filename);
    char *filename = strrchr(m->filename, '\\');
    if (filename == NULL) {
        filename = m->filename;
    }
    if (use_existing) {
        snprintf(objname, sizeof(objname), "%s%d", filename + 1, filenumber);
        fm_objname = objname;
    }
    DWORD alloc_size = getregionsize();
    DWORD alloc_mask = ~(alloc_size - 1);
    DWORD alloc_offset_mask = alloc_size - 1;
    if (filenumber == m->filenumber) { /* in the file denoted by m , must
                                          advance to page containing offset */
        m2->file_pos = offset;
        HANDLE fmhandle =
            CreateFileMapping(m->fh, NULL, PAGE_READWRITE, 0, 0, fm_objname);
        if (fmhandle == NULL) {
            return PH_ERRMMAPCREATE;
        }
        DWORD alloc_unit = alloc_mask & offset;
        DWORD alloc_offset = alloc_offset_mask & offset;
        DWORD offset_hi = 0, offset_lo = alloc_unit;
        m2->buf = (char *)MapViewOfFile(fmhandle, FILE_MAP_WRITE, offset_hi,
                                        offset_lo, alloc_offset + m->pgsize);
        if (m2->buf == NULL) {
            return PH_ERRMMAPVIEW;
        }
        m2->buf += alloc_offset;
        m2->hashdist = m->hashdist;
        m2->fh = m->fh;
        m2->branchfactor = m->branchfactor;
        m2->hash_type = m->hash_type;
        m2->pathlength = m->pathlength;
        m2->leafcapacity = m->leafcapacity;
        m2->nbdbfiles = m->nbdbfiles;
        m2->pgsize = m->pgsize;
        m2->filenumber = filenumber;
    } else { /* open and map to new file denoted by m->filename and filenumber
              */
        char extfile[256];
        snprintf(extfile, sizeof(extfile), "%s%d.mvp", m->filename, filenumber);
        m2->fh = CreateFile(extfile, GENERIC_READ | GENERIC_WRITE,
                            FILE_SHARE_READ | FILE_SHARE_WRITE, NULL,
                            OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
        if (m2->fh == INVALID_HANDLE_VALUE) {
            return PH_ERRFILEOPEN;
        }
        HANDLE fmhandle =
            CreateFileMapping(m2->fh, NULL, PAGE_READWRITE, 0, 0, fm_objname);
        if (fmhandle == NULL) {
            return PH_ERRMMAPCREATE;
        }
        m2->pathlength = m->pathlength;
        m2->hash_type = m->hash_type;
        m2->hashdist = m->hashdist;
        m2->branchfactor = m->branchfactor;
        m2->leafcapacity = m->leafcapacity;
        m2->nbdbfiles = m->nbdbfiles;
        m2->pgsize = m->pgsize;
        m2->file_pos = offset;
        m2->filenumber = filenumber;
        DWORD alloc_unit = alloc_mask & offset;
        DWORD alloc_offset = alloc_offset_mask & offset;
        DWORD offset_hi = 0, offset_lo = alloc_unit;
        m2->buf =
            (char *)MapViewOfFile(fmhandle, FILE_MAP_ALL_ACCESS, offset_hi,
                                  offset_lo, alloc_offset + m->pgsize);
        if (m2->buf == NULL) {
            return PH_ERRMMAPVIEW;
        }
        m2->buf += alloc_offset;
    }
    return PH_SUCCESS;
}

MVPRetCode _ph_unmap_mvpfile(uint8_t filenumber, off_t orig_pos, MVPFile *m,
                             MVPFile *m2) {
    sfree(m2->filename);
    if (filenumber == m->filenumber) { /*remap to same main file  */
        if (!UnmapViewOfFile(m2->buf)) {
            return PH_ERRMMAPUNVIEW;
        }
    } else { /*remap to back to main file  */
        if (!UnmapViewOfFile(m2->buf)) {
            return PH_ERRMMAPUNVIEW;
        }
        if (!CloseHandle(m2->fh)) {
            return PH_ERRFILECLOSE;
        }
    }
    return PH_SUCCESS;
}

__declspec(dllexport) float hammingdistance(DP *pntA, DP *pntB) {
    uint8_t htypeA = pntA->hash_type;
    uint8_t htypeB = pntB->hash_type;
    if (htypeA != htypeB) return -1.0;
    if (htypeA != UINT64ARRAY) return -1.0;
    if ((pntA->hash_length > 1) || (pntB->hash_length > 1)) return -1.0;
    ulong64 *hashA = (ulong64 *)pntA->hash;
    ulong64 *hashB = (ulong64 *)pntB->hash;
    int res = ph_hamming_distance(*hashA, *hashB);
    return (float)res;
}

__declspec(dllexport) MVPRetCode
    ph_query_mvptree(MVPFile *m, DP *query, int knearest, float radius,
                     float threshold, DP **results, int *count, int level) {
    int BranchFactor = m->branchfactor;
    int LengthM1 = BranchFactor - 1;
    int LengthM2 = BranchFactor * LengthM1;
    int PathLength = m->pathlength;
    MVPRetCode ret = PH_SUCCESS;

    if ((!m) || (!query) || (!results) || (!count) || (level < 0))
        return PH_ERRNULLARG;

    hash_compareCB hashdist = m->hashdist;

    if (!hashdist) return PH_ERRNULLARG;

    off_t offset_mask, page_mask;

    offset_mask = m->pgsize - 1;
    page_mask = ~(m->pgsize - 1);

    uint8_t ntype;
    memcpy(&ntype, &m->buf[m->file_pos & offset_mask], sizeof(uint8_t));
    m->file_pos++;
    if (ntype == 0) { /* leaf */
        DP *sv1 = ph_read_datapoint(m);
        DP *sv2 = ph_read_datapoint(m);
        if (sv1) {
            float d1 = hashdist(query, sv1);
            /* check if distance(sv1,query) <= radius  */
            if (d1 <= threshold) {
                results[(*count)++] = sv1;
                if (*count >= knearest) {
                    return PH_ERRCAP;
                }
            } else {
                hfree(sv1->id);
                hfree(sv1->hash);
                ph_free_datapoint(sv1);
            }
            if (sv2) {
                float d2 = hashdist(query, sv2);
                /* check if distance(sv2,query) <= radius */
                if (d2 <= threshold) {
                    results[(*count)++] = sv2;
                    if (*count >= knearest) {
                        return PH_ERRCAP;
                    }
                } else {
                    hfree(sv2->id);
                    hfree(sv2->hash);
                    ph_free_datapoint(sv2);
                }
                if (level < PathLength) {
                    query->path[level] = d1;
                }
                if (level + 1 < PathLength) {
                    query->path[level + 1] = d2;
                }
                uint8_t Np;
                memcpy(&Np, &m->buf[m->file_pos & offset_mask],
                       sizeof(uint8_t));
                m->file_pos += sizeof(uint8_t);
                off_t curr_pos;
                int include;
                /* read in each datapoint in the leaf - only retrieve the point
                if it dist(sv1,dp)=da  and dist(sv2,dp)=db cannot preclude the
                point */
                for (int i = 0; i < Np; i++) {
                    include = 1;
                    float da, db;
                    off_t point_offset;
                    memcpy(&da, &m->buf[m->file_pos & offset_mask],
                           sizeof(float));
                    m->file_pos += sizeof(float);
                    memcpy(&db, &m->buf[m->file_pos & offset_mask],
                           sizeof(float));
                    m->file_pos += sizeof(float);
                    if ((d1 - radius <= da) && (d1 + radius >= da) &&
                        (d2 - radius <= db) && (d2 + radius >= db)) {
                        // if (1) {
                        memcpy(&point_offset,
                               &m->buf[m->file_pos & offset_mask],
                               sizeof(off_t));
                        m->file_pos += sizeof(off_t);
                        curr_pos = m->file_pos;
                        m->file_pos = point_offset;
                        DP *dp = ph_read_datapoint(m);
                        m->file_pos = curr_pos;
                        /* test each path[] distance and as soon as one does not
                        fit disclude the point */
                        if (dp) {
                            int pl = (level <= PathLength) ? level : PathLength;
                            for (int j = 0; j < pl; j++) {
                                if (!((query->path[j] - radius <=
                                       dp->path[j]) &&
                                      (query->path[j] + radius >=
                                       dp->path[j]))) {
                                    include = 0;
                                    break;
                                }
                            }
                            if ((include) &&
                                (hashdist(query, dp) <= threshold)) {
                                results[(*count)++] = dp;
                            } else {
                                hfree(dp->id);
                                hfree(dp->hash);
                                ph_free_datapoint(dp);
                            }
                            if (*count >= knearest) {
                                return PH_ERRCAP;
                            }
                        }
                    } else {
                        m->file_pos += sizeof(off_t);
                    }
                }
            }
        }
    } else if (ntype == 1) { /* internal */
        /* read sv1, sv2 */
        DP *sv1 = ph_read_datapoint(m);
        DP *sv2 = ph_read_datapoint(m);

        /* read 1st and 2nd level pivots */
        float *M1 = (float *)malloc(LengthM1 * sizeof(float));
        if (!M1) return PH_ERRMEMALLOC;
        float *M2 = (float *)malloc(LengthM2 * sizeof(float));
        if (!M2) return PH_ERRMEMALLOC;

        memcpy(M1, &m->buf[m->file_pos & offset_mask],
               LengthM1 * sizeof(float));
        m->file_pos += LengthM1 * sizeof(float);
        memcpy(M2, &m->buf[m->file_pos & offset_mask],
               LengthM2 * sizeof(float));
        m->file_pos += LengthM2 * sizeof(float);

        float d1 = hashdist(query, sv1);
        float d2 = hashdist(query, sv2);

        /* fill in path values in query */
        if (level < PathLength) query->path[level] = d1;
        if (level + 1 < PathLength) query->path[level + 1] = d2;

        /* check if sv1 sv2 are close enough to query  */
        if (d1 <= threshold) {
            results[(*count)++] = sv1;
            if (*count >= knearest) {
                return PH_ERRCAP;
            }
        } else {
            hfree(sv1->id);
            hfree(sv1->hash);
            ph_free_datapoint(sv1);
        }

        if (d2 <= threshold) {
            results[(*count)++] = sv2;
            if (*count >= knearest) {
                return PH_ERRCAP;
            }
        } else {
            hfree(sv2->id);
            hfree(sv2->hash);
            ph_free_datapoint(sv2);
        }
        /* based on d1,d2 values, find appropriate child nodes to explore */

        int pivot1, pivot2;
        uint8_t filenumber;
        off_t child_pos, curr_pos;
        off_t start_pos = m->file_pos;
        /* check <= each M1 pivot */
        for (pivot1 = 0; pivot1 < LengthM1; pivot1++) {
            if (d1 - radius <= M1[pivot1]) {
                /* check <= each M2 pivot */
                for (pivot2 = 0; pivot2 < LengthM1; pivot2++) {
                    if (d2 - radius <= M2[pivot2 + pivot1 * LengthM1]) {
                        /*determine pos from which to read filenumber and offset
                         */
                        curr_pos =
                            start_pos + (pivot2 + pivot1 * BranchFactor) *
                                            (sizeof(uint8_t) + sizeof(off_t));
                        m->file_pos = curr_pos;
                        memcpy(&filenumber,
                               &(m->buf[m->file_pos & offset_mask]),
                               sizeof(uint8_t));
                        m->file_pos++;
                        memcpy(&child_pos, &(m->buf[m->file_pos & offset_mask]),
                               sizeof(off_t));
                        m->file_pos += sizeof(off_t);
                        /*save position and remap to new file/position  */
                        MVPFile m2;
                        if ((filenumber == 0) && (child_pos == 0)) {
                            continue;
                        }
                        ret = _ph_map_mvpfile(filenumber, child_pos, m, &m2, 1);
                        if (ret == PH_SUCCESS) {
                            ret = ph_query_mvptree(&m2, query, knearest, radius,
                                                   threshold, results, count,
                                                   level + 2);
                            if (ret != PH_SUCCESS) goto query_cleanup;
                            /* unmap and remap to the origional file/posion */
                            ret = _ph_unmap_mvpfile(filenumber, m->file_pos, m,
                                                    &m2);
                            if (ret != PH_SUCCESS) goto query_cleanup;
                        } else {
                            goto query_cleanup;
                        }
                    }
                }
                /* check > last M2 */
                if (d2 + radius >= M2[LengthM1 - 1 + pivot1 * LengthM1]) {
                    /*determine position from which to read filenumber and
                     * offset */
                    curr_pos =
                        start_pos + (BranchFactor - 1 + pivot1 * BranchFactor) *
                                        (sizeof(uint8_t) + sizeof(off_t));
                    m->file_pos = curr_pos;
                    memcpy(&filenumber, &m->buf[m->file_pos & offset_mask],
                           sizeof(uint8_t));
                    m->file_pos++;
                    memcpy(&child_pos, &m->buf[m->file_pos & offset_mask],
                           sizeof(off_t));
                    m->file_pos += sizeof(off_t);
                    /*saveposition and remap to new file/position */
                    MVPFile m2;
                    if ((filenumber == 0) && (child_pos == 0)) {
                        continue;
                    }
                    ret = _ph_map_mvpfile(filenumber, child_pos, m, &m2, 1);
                    if (ret == PH_SUCCESS) {
                        ret = ph_query_mvptree(&m2, query, knearest, radius,
                                               threshold, results, count,
                                               level + 2);
                        if (ret != PH_SUCCESS) goto query_cleanup;
                        /*unmap and remap to original file/position  */
                        ret =
                            _ph_unmap_mvpfile(filenumber, m->file_pos, m, &m2);
                        if (ret != PH_SUCCESS) goto query_cleanup;
                    } else {
                        goto query_cleanup;
                    }
                }
            }
        }
        /* check >=  last M1 pivot */
        if (d1 + radius >= M1[LengthM1 - 1]) {
            /* check <= each M2 pivot */
            for (pivot2 = 0; pivot2 < LengthM1; pivot2++) {
                if (d2 - radius <= M2[pivot2 + LengthM1 * LengthM1]) {
                    /*determine pos from which to read filenumber and position
                     */
                    curr_pos =
                        start_pos + (pivot2 + LengthM1 * BranchFactor) *
                                        (sizeof(uint8_t) + sizeof(off_t));
                    m->file_pos = curr_pos;
                    memcpy(&filenumber, &m->buf[m->file_pos & offset_mask],
                           sizeof(uint8_t));
                    m->file_pos++;
                    memcpy(&child_pos, &m->buf[m->file_pos & offset_mask],
                           sizeof(off_t));
                    m->file_pos += sizeof(off_t);
                    /*save file position and remap to new filenumber/offset  */
                    MVPFile m2;
                    if ((filenumber == 0) && (child_pos == 0)) {
                        continue;
                    }
                    ret = _ph_map_mvpfile(filenumber, child_pos, m, &m2, 1);
                    if (ret == PH_SUCCESS) {
                        ret = ph_query_mvptree(&m2, query, knearest, radius,
                                               threshold, results, count,
                                               level + 2);
                        if (ret != PH_SUCCESS) goto query_cleanup;
                        /* unmap/remap to original filenumber/position */
                        ret =
                            _ph_unmap_mvpfile(filenumber, m->file_pos, m, &m2);
                        if (ret != PH_SUCCESS) goto query_cleanup;
                    } else {
                        goto query_cleanup;
                    }
                }
            }

            /* check >= last M2 pivot */
            if (d2 + radius >= M2[LengthM1 - 1 + LengthM1 * LengthM1]) {
                /* determine position from which to read filenumber and child
                 * position */
                curr_pos =
                    start_pos + (BranchFactor - 1 + LengthM1 * BranchFactor) *
                                    (sizeof(uint8_t) + sizeof(off_t));
                m->file_pos = curr_pos;
                memcpy(&filenumber, &m->buf[m->file_pos & offset_mask],
                       sizeof(uint8_t));
                m->file_pos += sizeof(uint8_t);
                memcpy(&child_pos, &m->buf[m->file_pos & offset_mask],
                       sizeof(off_t));
                m->file_pos += sizeof(off_t);
                /* save position and remap to new filenumber/position */
                MVPFile m2;
                if (!(filenumber == 0 && child_pos == 0)) {
                    ret = _ph_map_mvpfile(filenumber, child_pos, m, &m2, 1);
                    if (ret == PH_SUCCESS) {
                        ret = ph_query_mvptree(&m2, query, knearest, radius,
                                               threshold, results, count,
                                               level + 2);
                        if (ret != PH_SUCCESS) goto query_cleanup;
                        /* return to original and remap to original
                         * filenumber/position */
                        ret =
                            _ph_unmap_mvpfile(filenumber, m->file_pos, m, &m2);
                        if (ret != PH_SUCCESS) goto query_cleanup;
                    } else {
                        goto query_cleanup;
                    }
                }
            }
        }
    query_cleanup:
        hfree(M1);
        hfree(M2);
    } else { /* unrecognized node */
        ret = PH_ERRNTYPE;
    }
    return ret;
}

__declspec(dllexport) MVPRetCode
    ph_query_mvptree(MVPFile *m, DP *query, int knearest, float radius,
                     float threshold, DP **results, int *count) {
    /*use host pg size until file pg size used can be determined  */
    m->pgsize = (off_t)getpagesize();

    char mainfile[MAX_PATH];
    snprintf(mainfile, sizeof(mainfile), "%s.mvp", m->filename);
    m->fh = CreateFile(mainfile, GENERIC_READ | GENERIC_WRITE,
                       FILE_SHARE_READ | FILE_SHARE_WRITE, NULL, OPEN_EXISTING,
                       FILE_ATTRIBUTE_NORMAL, NULL);
    if (m->fh == INVALID_HANDLE_VALUE) {
        return PH_ERRFILEOPEN;
    }
    char *filename = strrchr(m->filename, '\\');
    if (filename == NULL) {
        filename = m->filename;
    }
    char fm_objname[32];
    snprintf(fm_objname, sizeof(fm_objname), "%s%d", filename + 1, 0);

    m->file_pos = 0;
    HANDLE fmhandle =
        CreateFileMapping(m->fh, NULL, PAGE_READWRITE, 0, 0, fm_objname);
    if (fmhandle == NULL) {
        return PH_ERRMMAPCREATE;
    }
    m->buf =
        (char *)MapViewOfFile(fmhandle, FILE_MAP_ALL_ACCESS, 0, 0, m->pgsize);
    if (m->buf == NULL) {
        return PH_ERRMMAPVIEW;
    }
    char tag[17];
    int version;
    int int_pgsize, leaf_pgsize;
    uint8_t nbdbfiles, bf, p, k, type;

    memcpy((char *)tag, (char *)&(m->buf[m->file_pos]), 16);
    tag[16] = '\0';
    m->file_pos += 16;

    memcpy(&version, &m->buf[m->file_pos], sizeof(int));
    m->file_pos += sizeof(int);

    memcpy(&int_pgsize, &m->buf[m->file_pos], sizeof(int));
    m->file_pos += sizeof(int);

    memcpy(&leaf_pgsize, &m->buf[m->file_pos], sizeof(int));
    m->file_pos += sizeof(int);

    memcpy(&nbdbfiles, &m->buf[m->file_pos++], 1);

    memcpy(&bf, &m->buf[m->file_pos++], 1);

    memcpy(&p, &m->buf[m->file_pos++], 1);

    memcpy(&k, &m->buf[m->file_pos++], 1);

    memcpy(&type, &m->buf[m->file_pos++], 1);

    if (!UnmapViewOfFile(m->buf)) {
        return PH_ERRMMAPUNVIEW;
    }
    m->buf =
        (char *)MapViewOfFile(fmhandle, FILE_MAP_ALL_ACCESS, 0, 0, int_pgsize);
    if (m->buf == NULL) {
        return PH_ERRMMAPVIEW;
    }
    m->branchfactor = bf;
    m->pathlength = p;
    m->leafcapacity = k;
    m->pgsize = (off_t)int_pgsize;
    m->hash_type = (HashType)type;
    m->file_pos = HeaderSize;
    m->nbdbfiles = nbdbfiles;
    m->filenumber = 0;

    /* add path onto the query */
    if (query->path == NULL) {
        query->path = (float *)malloc((int)(p) * sizeof(float));
    }

    /* finish the query by calling the recursive auxiliary function */
    *count = 0;
    MVPRetCode res = ph_query_mvptree(m, query, knearest, radius, threshold,
                                      results, count, 0);

    if (!UnmapViewOfFile(m->buf)) {
        res = PH_ERRMMAPUNVIEW;
    }
    if (!CloseHandle(fmhandle)) {
        res = PH_ERRFILECLOSE;
    }
    if (!CloseHandle(m->fh)) {
        res = PH_ERRFILECLOSE;
    }
    hfree(query->path);

    for (int i = 0; i < *count; i++) {
        hfree(results[i]->path);
    }

    return res;
}

__declspec(dllexport) MVPRetCode
    ph_save_mvptree(MVPFile *m, DP **points, int nbpoints, int saveall_flag,
                    int level, FileIndex *pOffset) {
    int Np = (nbpoints >= 2) ? nbpoints - 2 : 0;
    int BranchFactor = m->branchfactor;
    int PathLength = m->pathlength;
    int LeafCapacity = m->leafcapacity;
    int LengthM1 = BranchFactor - 1;
    int LengthM2 = BranchFactor * LengthM1;
    int Fanout = BranchFactor * BranchFactor;
    if ((!m) || (!points) || (nbpoints < 0) || (pOffset == NULL)) {
        return PH_ERRNULLARG;
    }
    MVPRetCode ret = PH_SUCCESS;
    hash_compareCB hashdist = m->hashdist;
    if (hashdist == NULL) {
        return PH_ERRNULLARG;
    }
    DWORD alloc_size = getregionsize();
    DWORD alloc_mask = ~(alloc_size - 1);
    DWORD alloc_offset_mask = alloc_size - 1;
    off_t offset_mask = m->pgsize - 1;
    off_t page_mask = ~(m->pgsize - 1);
    if (nbpoints <= 0) {
        pOffset->fileno = 0;
        pOffset->offset = 0;
        return PH_SUCCESS;
    } else if (nbpoints <= LeafCapacity + 2) { /* leaf */
        uint8_t ntype = 0;
        MVPFile m2;
        m2.branchfactor = m->branchfactor;
        m2.leafcapacity = m->leafcapacity;
        m2.pathlength = m->pathlength;
        m2.pgsize = m->pgsize;
        m2.hashdist = hashdist;
        m2.hash_type = m->hash_type;
        /* open new file */
        char extfile[256];
        snprintf(extfile, sizeof(extfile), "%s%d.mvp", m->filename,
                 m->nbdbfiles);
        m2.fh =
            CreateFile(extfile, GENERIC_READ | GENERIC_WRITE,
                       FILE_SHARE_WRITE | FILE_SHARE_READ | FILE_SHARE_DELETE,
                       NULL, OPEN_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
        if (m2.fh == INVALID_HANDLE_VALUE) {
            free(pOffset);
            return PH_ERRFILEOPEN;
        }
        DWORD fsize = GetFileSize(m2.fh, NULL);
        /* test file size and open new file if necessary */
        while (fsize >= MaxFileSize) {
            CloseHandle(m2.fh);
            m->nbdbfiles++;
            snprintf(extfile, sizeof(extfile), "%s%d.mvp", m->filename,
                     m->nbdbfiles);
            m2.fh = CreateFile(
                extfile, GENERIC_READ | GENERIC_WRITE,
                FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE, NULL,
                OPEN_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
            if (m2.fh < INVALID_HANDLE_VALUE) {
                free(pOffset);
                return PH_ERRFILEOPEN;
            }
            fsize = GetFileSize(m2.fh, NULL);
        }

        char fm_objname[256];
        snprintf(fm_objname, sizeof(fm_objname), "%s%d", m->filename,
                 m->nbdbfiles);

        m2.file_pos = fsize;
        off_t end_pos = fsize + m->pgsize;
        pOffset->fileno = m->nbdbfiles;
        pOffset->offset = m2.file_pos;
        off_t pa_offset = m2.file_pos & page_mask;

        if (SetFilePointer(m2.fh, m->pgsize, 0, FILE_END) ==
            INVALID_SET_FILE_POINTER) {
            return PH_ERRFILESETPOINTER;
        }
        if (!SetEndOfFile(m2.fh)) {
            return PH_ERRFILESETEOF;
        }

        HANDLE fmhandle = CreateFileMapping(m2.fh, NULL, PAGE_READWRITE, 0,
                                            (DWORD)end_pos, NULL);
        if (fmhandle == NULL) {
            return PH_ERRMMAPCREATE;
        }

        DWORD alloc_unit = alloc_mask & m2.file_pos;
        DWORD alloc_offset = alloc_offset_mask & m2.file_pos;

        DWORD offset_low = alloc_unit;
        DWORD offset_hi = 0;

        m2.buf = (char *)MapViewOfFile(fmhandle, FILE_MAP_ALL_ACCESS, offset_hi,
                                       offset_low, alloc_offset + m2.pgsize);
        if (m2.buf == NULL) {
            return PH_ERRMMAPVIEW;
        }
        m2.buf += alloc_offset;

        m2.buf[m2.file_pos & offset_mask] = ntype;
        m2.file_pos++;
        /* find vantage points, sv1 and sv2 */
        int sv1_pos, sv2_pos;
        ph_selectvantagepoints(&m2, points, nbpoints, sv1_pos, sv2_pos);
        DP *sv1 = NULL, *sv2 = NULL;
        if (sv1_pos >= 0) {
            sv1 = points[sv1_pos];
        }
        if (sv2_pos >= 0) {
            sv2 = points[sv2_pos];
        }
        /* if file pos is beyond pg size*/
        if ((m2.file_pos & offset_mask) + ph_sizeof_dp(sv1, &m2) > m->pgsize)
            return PH_SMPGSIZE;

        ph_save_datapoint(sv1, &m2);

        if (!sv1) goto leafcleanup;

        /*if file pos is beyond pg size */
        if ((m2.file_pos & offset_mask) + ph_sizeof_dp(sv2, &m2) > m->pgsize) {
            return PH_SMPGSIZE;
        }

        ph_save_datapoint(sv2, &m2);

        if (!sv2) goto leafcleanup;

        m2.buf[m2.file_pos & offset_mask] = (uint8_t)Np;
        m2.file_pos++;

        /* write the leaf points */
        float d1, d2;
        off_t curr_pos = m2.file_pos;
        off_t last_pos =
            curr_pos + LeafCapacity * (2 * sizeof(float) + sizeof(off_t));
        off_t dp_pos;
        for (int i = 0; i < nbpoints; i++) {
            if ((i == sv1_pos) || (i == sv2_pos)) /* skip sv1 and sv2 */
                continue;
            /* write d1, d2 */
            d1 = hashdist(sv1, points[i]);
            d2 = hashdist(sv2, points[i]);
            if (level < PathLength) {
                points[i]->path[level] = d1;
            }
            if (level + 1 < PathLength) {
                points[i]->path[level + 1] = d2;
            }
            memcpy(&(m2.buf[m2.file_pos & offset_mask]), &d1, sizeof(float));
            m2.file_pos += sizeof(float);
            memcpy(&(m2.buf[m2.file_pos & offset_mask]), &d2, sizeof(float));
            m2.file_pos += sizeof(float);

            /* write the point[i] at the end and return to write what the offset
             * is*/
            curr_pos = m2.file_pos;
            m2.file_pos = last_pos;

            /* if file pos is beyond pg size */
            if ((m2.file_pos & offset_mask) + ph_sizeof_dp(points[i], &m2) >
                m->pgsize)
                return PH_SMPGSIZE;

            dp_pos = ph_save_datapoint(points[i], &m2);
            last_pos = m2.file_pos;
            m2.file_pos = curr_pos;

            memcpy(&(m2.buf[m2.file_pos & offset_mask]), &dp_pos,
                   sizeof(off_t));
            m2.file_pos += sizeof(off_t);
        }
    leafcleanup:
        FlushViewOfFile(m2.buf, m2.pgsize);
        UnmapViewOfFile(m2.buf);
        CloseHandle(fmhandle);
        CloseHandle(m2.fh);
    } else { /* internal node */
        pOffset->fileno = 0;
        pOffset->offset = m->file_pos;
        off_t orig_pos = m->file_pos;
        off_t end_pos = orig_pos + m->pgsize;
        char fm_objname[256];
        snprintf(fm_objname, sizeof(fm_objname), "%s%d", m->filename, 0);
        HANDLE fmhandle;
        char *buf = m->buf;
        if ((level > 0) &&
            (saveall_flag ==
             1)) { /* append new page to mainfile, unmap/map to it */
            DWORD fsize = GetFileSize(m->fh, NULL);
            m->file_pos = fsize;
            end_pos = m->file_pos + m->pgsize;

            if (SetFilePointer(m->fh, m->pgsize, 0, FILE_END) ==
                INVALID_SET_FILE_POINTER) {
                return PH_ERRFILESETPOINTER;
            }
            if (!SetEndOfFile(m->fh)) {
                return PH_ERRFILESETEOF;
            }
            fmhandle =
                CreateFileMapping(m->fh, NULL, PAGE_READWRITE, 0, 0, NULL);
            if (fmhandle == NULL) {
                return PH_ERRMMAPCREATE;
            }
            DWORD alloc_unit = alloc_mask & m->file_pos;
            DWORD alloc_offset = alloc_offset_mask & m->file_pos;
            DWORD offset_hi = 0, offset_lo = alloc_unit;
            m->buf =
                (char *)MapViewOfFile(fmhandle, FILE_MAP_ALL_ACCESS, offset_hi,
                                      offset_lo, alloc_offset + m->pgsize);
            if (m->buf == NULL) {
                return PH_ERRMMAPVIEW;
            }
            m->buf += alloc_offset;
            pOffset->offset = m->file_pos;
        }

        uint8_t ntype = 1;
        memcpy(&m->buf[m->file_pos++ & offset_mask], &ntype, sizeof(uint8_t));

        /* choose vantage points, sv1, sv2 */
        int sv1_pos, sv2_pos;
        float max_distance = 0.0f, min_distance = (float)INT_MAX;
        ph_selectvantagepoints(m, points, nbpoints, sv1_pos, sv2_pos);
        DP *sv1 = points[sv1_pos];
        DP *sv2 = points[sv2_pos];

        /* save sv1, sv2 */
        ph_save_datapoint(sv1, m);
        ph_save_datapoint(sv2, m);

        for (int i = 0; i < nbpoints; i++) {
            if (i != sv1_pos) {
                float d = m->hashdist(points[i], points[sv1_pos]);
                if (d > max_distance) {
                    max_distance = d;
                }
                if (d < min_distance) {
                    min_distance = d;
                }
            }
        }

        /* 1st tier pivots, M1, derived from the distance of each point from
         * sv1*/
        float step = (max_distance - min_distance) / BranchFactor;
        if (step <= 0.001) {
            return PH_ERRDISTFUNC;
        }

        float incr = step;

        float *M1 = (float *)malloc(LengthM1 * sizeof(float));
        float *M2 = (float *)malloc(LengthM2 * sizeof(float));
        if (!M1 || !M2) {
            return PH_ERRMEMALLOC;
        }

        for (int i = 0; i < LengthM1; i++) {
            M1[i] = min_distance + incr;
            incr += step;
            memcpy(&(m->buf[m->file_pos & offset_mask]), &M1[i], sizeof(float));
            m->file_pos += sizeof(float);
        }

        /*set up 1st tier sorting bins - contain pointers to DP points[] param
  so that we only move pointers, not the actual datapoints */
        DP ***bins = (DP ***)malloc(BranchFactor * sizeof(DP **));
        if (!bins) {
            hfree(M1);
            hfree(M2);
            return PH_ERRMEMALLOC;
        }
        int *mlens = (int *)calloc(BranchFactor,
                                   sizeof(int)); /*no. points in each bin */
        if (!mlens) {
            hfree(M1);
            hfree(M2);
            hfree(bins);
            return PH_ERRMEMALLOC;
        }

        for (int i = 0; i < BranchFactor; i++) {
            bins[i] = (DP **)malloc(
                Np * sizeof(DP **)); /*Np should be more than enough */
            if (!bins[i]) {
                hfree(M1);
                hfree(M2);
                hfree(bins);
                cfree(mlens);
                return PH_ERRMEMALLOC;
            }
        }

        /* sort points into bins (except sv1 and sv2 )*/
        for (int i = 0; i < nbpoints; i++) {
            if ((i == sv1_pos) || (i == sv2_pos)) continue;
            float cur_dist = m->hashdist(sv1, points[i]);
            if (level < PathLength) {
                points[i]->path[level] = cur_dist;
            }
            /* check if <= M1[i] */
            for (int j = 0; j < LengthM1; j++) {
                if (cur_dist <= M1[j]) {
                    bins[j][mlens[j]] = points[i];
                    mlens[j]++;
                    break;
                }
            }
            /* check if > last M1[] pivot */
            if (cur_dist > M1[LengthM1 - 1]) {
                bins[BranchFactor - 1][mlens[BranchFactor - 1]] = points[i];
                mlens[BranchFactor - 1]++;
            }
        }

        /* print 1st level sort bins
        for (int i=0; i<BranchFactor;i++){
            int row_len = mlens[i];
                for (int j=0;j<row_len;j++){
                    printf(" %d %d %s\n", i, j, bins[i][j]->id);
                }
        }
        */

        /* set up 2nd tier sorting bins */
        /* each row from bins to be sorted into bins2, each in turn */
        DP ***bins2 = (DP ***)malloc(BranchFactor * sizeof(DP ***));
        if (!bins2) {
            hfree(M1);
            hfree(M2);
            hfree(bins);
            cfree(mlens);
            return PH_ERRMEMALLOC;
        }
        int *mlens2 = (int *)calloc(BranchFactor,
                                    sizeof(int)); /*number points in each bin */
        if (!mlens2) {
            hfree(M1);
            hfree(M2);
            hfree(bins);
            hfree(bins2);
            cfree(mlens);
            return PH_ERRMEMALLOC;
        }

        for (int i = 0; i < BranchFactor; i++) {
            bins2[i] =
                (DP **)malloc(Np * sizeof(DP **)); /* Np is more than enough */
            if (!bins2[i]) {
                hfree(M1);
                hfree(M2);
                hfree(bins);
                cfree(mlens);
                hfree(bins2);
                cfree(mlens2);
                return PH_ERRMEMALLOC;
            }
        }
        FileIndex ChildIndex;
        off_t m2_pos = m->file_pos; /* start of M2 pivots */
        off_t child_pos =
            m2_pos +
            LengthM2 * sizeof(float); /*pos where child offsets are written*/
        off_t last_pos =
            child_pos + Fanout * (sizeof(uint8_t) + sizeof(off_t)); /* last pos
                                                   in in internal node */
        float *distance_vector = NULL;
        /* for each row of bin, sort the row into bins2 */
        for (int i = 0; i < BranchFactor; i++) {
            int row_len = mlens[i]; /* length of current row, bins[i] */
            for (int j = 0; j < BranchFactor;
                 j++) { /* reset the lengths to 0 */
                mlens2[j] = 0;
            }
            float *distance_vector = (float *)malloc(row_len * sizeof(float));
            if (!distance_vector) {
                hfree(M1);
                hfree(M2);
                hfree(bins);
                cfree(mlens);
                hfree(bins2);
                cfree(mlens2);
                return PH_ERRMEMALLOC;
            }

            /* 2nd tier pivots M2[], for row */
            max_distance = 0;
            min_distance = 100000000.0f;
            for (int j = 0; j < row_len; j++) {
                distance_vector[j] = hashdist(sv2, bins[i][j]);
                if (distance_vector[j] > max_distance)
                    max_distance = distance_vector[j];
                if (distance_vector[j] < min_distance)
                    min_distance = distance_vector[j];
                if (level + 1 < PathLength) {
                    bins[i][j]->path[level + 1] = distance_vector[j];
                }
            }

            step = (max_distance - min_distance) / BranchFactor;
            incr = step;

            for (int j = 0; j < LengthM1; j++) {
                M2[j + i * LengthM1] = min_distance + incr;
                memcpy(&(m->buf[m2_pos & offset_mask]), &M2[j + i * LengthM1],
                       sizeof(float));
                incr += step;
                m2_pos += sizeof(float);
            }

            /* sort bins[i] into bins2 */
            for (int j = 0; j < row_len; j++) {
                DP *current = bins[i][j];
                /*check <= each M2 pivot  */
                for (int k = 0; k < LengthM1; k++) {
                    if (distance_vector[j] <= M2[k + i * LengthM1]) {
                        bins2[k][mlens2[k]] = current;
                        mlens2[k]++;
                        break;
                    }
                }
                /* check > last M2 pivot  */
                if (distance_vector[j] > M2[LengthM1 - 1 + i * LengthM1]) {
                    bins2[BranchFactor - 1][mlens2[BranchFactor - 1]] = current;
                    mlens2[BranchFactor - 1]++;
                }
            }

            /* print 2nd tier sort bins
            for (int j=0;j<BranchFactor;j++){
                int r2 = mlens2[j];
                for (int k=0;k<r2;k++){
                printf(" %d %d %d %s\n", i, j, k, bins2[j][k]->id);
                    }
            }
            */
            /* save child nodes */
            for (int j = 0; j < BranchFactor; j++) {
                m->file_pos = last_pos;
                ret = ph_save_mvptree(m, bins2[j], mlens2[j], saveall_flag,
                                      level + 2, &ChildIndex);
                if (ret == PH_SUCCESS) {
                    last_pos = m->file_pos;
                    m->file_pos = child_pos;
                    memcpy(&(m->buf[m->file_pos & offset_mask]),
                           &(ChildIndex.fileno), sizeof(uint8_t));
                    m->file_pos++;
                    memcpy(&(m->buf[m->file_pos & offset_mask]),
                           &(ChildIndex.offset), sizeof(off_t));
                    m->file_pos += sizeof(off_t);
                    child_pos = m->file_pos;
                } else {
                    goto cleanup;
                }
            }
            m->file_pos = last_pos;
            hfree(distance_vector);
        }

        /* remap to orig_pos */
        if ((level > 0) && (saveall_flag == 1)) {
            /* unmap/remap to page with original position */
            if (!FlushViewOfFile(m->buf, m->pgsize)) {
                return PH_ERRMMAPFLUSH;
            }
            if (!UnmapViewOfFile(m->buf)) {
                return PH_ERRMMAPUNVIEW;
            }
            if (!CloseHandle(fmhandle)) {
                return PH_ERRFILECLOSE;
            }
            m->buf = buf;
            m->file_pos = orig_pos;
        }

        /* cleanup */
    cleanup:
        hfree(bins);
        hfree(bins2);
        cfree(mlens);
        cfree(mlens2);
        hfree(distance_vector);
        hfree(M1);
        hfree(M2);
    }
    return ret;
}

__declspec(dllexport) MVPRetCode
    ph_save_mvptree(MVPFile *m, DP **points, int nbpoints) {
    if ((!m) || (!points)) {
        return PH_ERRARG;
    }

    if (m->pgsize == 0) /*use host pg size as default */
        m->pgsize = (off_t)getpagesize();

    /* check to see that the pg sizes are at least the size of host page size */
    off_t host_pgsize = (off_t)getpagesize();
    if (m->pgsize < host_pgsize) {
        return PH_ERRPGSZINCOMP;
    }

    /* pg sizes must be a power of 2 */
    if ((m->pgsize) & (m->pgsize - 1)) {
        return PH_ERRPGSZINCOMP;
    }
    if ((nbpoints < m->leafcapacity + 2) || (!m) || (!points)) {
        return PH_INSUFFPOINTS;
    }
    /* add paths to datapoints */
    for (int i = 0; i < nbpoints; i++) {
        points[i]->path = (float *)malloc((m->pathlength) * sizeof(float));
    }
    /* open main file */
    char mainfile[MAX_PATH];
    snprintf(mainfile, sizeof(mainfile), "%s.mvp", m->filename);
    m->fh = CreateFile(mainfile, GENERIC_READ | GENERIC_WRITE,
                       FILE_SHARE_WRITE | FILE_SHARE_READ | FILE_SHARE_DELETE,
                       NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
    if (m->fh == INVALID_HANDLE_VALUE) {
        return PH_ERRFILEOPEN;
    }
    /* enlarge to page*/
    m->file_pos = 0;
    if (SetFilePointer(m->fh, m->pgsize, 0, FILE_END) ==
        INVALID_SET_FILE_POINTER) {
        return PH_ERRFILESETPOINTER;
    }
    if (!SetEndOfFile(m->fh)) {
        return PH_ERRFILESETEOF;
    }
    char *filename = strrchr(m->filename, '\\');
    if (filename == NULL) {
        filename = m->filename;
    }
    char fm_objname[32];
    snprintf(fm_objname, sizeof(fm_objname), "%s%d", filename, 0);

    HANDLE fmhandle =
        CreateFileMapping(m->fh, NULL, PAGE_READWRITE, 0, 0, fm_objname);
    if (fmhandle == NULL) {
        return PH_ERRMMAPCREATE;
    }
    DWORD alloc_unit = 0;
    DWORD alloc_offset = 0;
    DWORD offset_hi = 0, offset_lo = alloc_unit;
    m->buf = (char *)MapViewOfFile(fmhandle, FILE_MAP_ALL_ACCESS, offset_hi,
                                   offset_lo, alloc_offset + m->pgsize);
    m->buf += alloc_offset;
    if (m->buf == NULL) {
        return PH_ERRMMAPVIEW;
    }

    m->nbdbfiles = 1;
    /* write header within first HeaderSize bytes */
    int version = 0;
    int int_pgsize = (int)(m->pgsize);
    int leaf_pgsize = (int)(m->pgsize);

    memcpy(&m->buf[m->file_pos], mvptag, 16);
    m->file_pos += 16;

    memcpy(&m->buf[m->file_pos], &version, sizeof(int));
    m->file_pos += sizeof(int);

    memcpy(&m->buf[m->file_pos], &int_pgsize, sizeof(int));
    m->file_pos += sizeof(int);

    memcpy(&m->buf[m->file_pos], &leaf_pgsize, sizeof(int));
    m->file_pos += sizeof(int);

    off_t nbfiles_pos = m->file_pos;
    memcpy(&m->buf[m->file_pos++], &m->nbdbfiles, sizeof(uint8_t));

    memcpy(&m->buf[m->file_pos++], &m->branchfactor, sizeof(uint8_t));

    memcpy(&m->buf[m->file_pos++], &m->pathlength, sizeof(uint8_t));

    memcpy(&m->buf[m->file_pos++], &m->leafcapacity, sizeof(uint8_t));
    uint8_t type = (uint8_t)m->hash_type;
    memcpy(&m->buf[m->file_pos++], &type, sizeof(uint8_t));

    m->file_pos = HeaderSize;

    FileIndex node_offset;
    MVPRetCode ret = PH_SUCCESS;
    ret = ph_save_mvptree(m, points, nbpoints, 1, 0, &node_offset);

    memcpy(&m->buf[nbfiles_pos], &m->nbdbfiles, sizeof(uint8_t));

    if (!FlushViewOfFile(m->buf, m->pgsize)) {
        return PH_ERRMMAPFLUSH;
    }

    if (!UnmapViewOfFile(m->buf)) {
        ret = PH_ERRMMAPUNVIEW;
        return ret;
    }
    if (!CloseHandle(fmhandle)) {
        ret = PH_ERRFILECLOSE;
    }
    if (!CloseHandle(m->fh)) {
        ret = PH_ERRFILECLOSE;
    }
    return ret;
}

__declspec(dllexport) MVPRetCode
    ph_add_mvptree(MVPFile *m, DP *new_dp, int level) {
    MVPRetCode ret = PH_SUCCESS;
    uint8_t ntype;
    off_t offset_mask, page_mask;

    if ((!m) || (!new_dp) || (level < 0)) return PH_ERRNULLARG;

    hash_compareCB hashdist = m->hashdist;

    if (!hashdist) return PH_ERRNULLARG;

    offset_mask = m->pgsize - 1;
    page_mask = ~(m->pgsize - 1);
    int LengthM1 = m->branchfactor - 1;
    int LengthM2 = m->branchfactor * LengthM1;

    off_t start_pos = m->file_pos;
    memcpy(&ntype, &m->buf[m->file_pos & offset_mask], sizeof(uint8_t));
    m->file_pos++;

    if (ntype == 0) {
        uint8_t Np = 0;
        DP *sv1 = ph_read_datapoint(m);
        if (sv1) {
            DP *sv2 = ph_read_datapoint(m);
            float d1 = hashdist(sv1, new_dp);
            if (level < m->pathlength) {
                new_dp->path[level] = d1;
            }
            if (sv2) {
                off_t Np_pos = m->file_pos;
                memcpy(&Np, &m->buf[m->file_pos & offset_mask],
                       sizeof(uint8_t));
                m->file_pos++;

                off_t offset_start = m->file_pos;
                float d2 = hashdist(sv2, new_dp);
                if (level + 1 < m->pathlength) {
                    new_dp->path[level + 1] = d2;
                }
                off_t curr_pos, new_pos, point_pos;
                if (Np == 0) { /* add new dp to existing leaf */
                    memcpy(&m->buf[m->file_pos & offset_mask], &d1,
                           sizeof(float));
                    m->file_pos += sizeof(float);
                    memcpy(&m->buf[m->file_pos & offset_mask], &d2,
                           sizeof(float));
                    m->file_pos += sizeof(float);

                    curr_pos = m->file_pos;

                    m->file_pos =
                        offset_start +
                        (m->leafcapacity) * (2 * sizeof(float) + sizeof(off_t));

                    if ((m->file_pos & offset_mask) + ph_sizeof_dp(new_dp, m) >
                        m->pgsize)
                        return PH_SMPGSIZE;

                    new_pos = ph_save_datapoint(new_dp, m);

                    m->file_pos = curr_pos;
                    memcpy(&m->buf[m->file_pos & offset_mask], &new_pos,
                           sizeof(off_t));

                    Np++;
                    memcpy(&m->buf[Np_pos & offset_mask], &Np, sizeof(uint8_t));

                } else if (Np <
                           m->leafcapacity) { /* add new dp to existing link */
                    m->file_pos +=
                        (Np - 1) * (2 * sizeof(float) + sizeof(off_t));
                    m->file_pos += 2 * sizeof(float);
                    memcpy(&point_pos, &m->buf[m->file_pos & offset_mask],
                           sizeof(off_t));
                    m->file_pos += sizeof(off_t);

                    curr_pos = m->file_pos;
                    m->file_pos = point_pos;
                    uint8_t active;
                    uint16_t byte_len;
                    memcpy(&active, &m->buf[m->file_pos & offset_mask],
                           sizeof(uint8_t));
                    m->file_pos++;
                    memcpy(&byte_len, &m->buf[m->file_pos & offset_mask],
                           sizeof(uint16_t));
                    m->file_pos += sizeof(uint16_t);
                    m->file_pos += byte_len;

                    if ((m->file_pos & offset_mask) + ph_sizeof_dp(new_dp, m) >
                        m->pgsize)
                        return PH_SMPGSIZE;

                    new_pos = ph_save_datapoint(new_dp, m);
                    memcpy(&m->buf[curr_pos & offset_mask], &d1, sizeof(float));
                    curr_pos += sizeof(float);
                    memcpy(&m->buf[curr_pos & offset_mask], &d2, sizeof(float));
                    curr_pos += sizeof(float);
                    memcpy(&m->buf[curr_pos & offset_mask], &new_pos,
                           sizeof(off_t));
                    curr_pos += sizeof(off_t);
                    Np++;
                    memcpy(&m->buf[Np_pos & offset_mask], &Np, sizeof(uint8_t));
                } else { /* convert to internal node and add leaf nodes */
                    DP **points = (DP **)malloc((Np + 3) * sizeof(DP **));
                    points[0] = sv1;
                    points[1] = sv2;

                    for (int i = 2; i < Np + 2; i++) {
                        m->file_pos += 2 * sizeof(float);
                        memcpy(&point_pos, &m->buf[m->file_pos & offset_mask],
                               sizeof(off_t));
                        m->file_pos += sizeof(off_t);

                        curr_pos = m->file_pos;
                        m->file_pos = point_pos;
                        points[i] = ph_read_datapoint(m);
                        m->file_pos = curr_pos;
                    }
                    points[Np + 2] = new_dp;
                    m->file_pos = start_pos;
                    FileIndex ChildOffset;
                    ret = ph_save_mvptree(m, points, Np + 3, 0, level,
                                          &ChildOffset);
                    if (ret != PH_SUCCESS) {
                        free(sv1->id);
                        free(sv1->hash);
                        ph_free_datapoint(sv1);
                        free(sv2->id);
                        free(sv2->hash);
                        ph_free_datapoint(sv2);
                        for (int i = 2; i < Np + 2; i++) {
                            hfree(points[i]->id);
                            hfree(points[i]->hash);
                            hfree(points[i]);
                        }
                        return ret;
                    }
                    for (int i = 2; i < Np + 2; i++) {
                        hfree(points[i]->id);
                        hfree(points[i]->hash);
                        ph_free_datapoint(points[i]);
                    }
                    hfree(points);
                }
                hfree(sv2->id);
                hfree(sv2->hash);
                ph_free_datapoint(sv2);
            } else { /* put new point into sv2 pos */
                m->file_pos = start_pos;
                ntype = 0;
                memcpy(&m->buf[m->file_pos & offset_mask], &ntype,
                       sizeof(uint8_t));
                m->file_pos++;

                if ((m->file_pos & offset_mask) + ph_sizeof_dp(sv1, m) >
                    m->pgsize)
                    return PH_SMPGSIZE;
                ph_save_datapoint(sv1, m);

                if ((m->file_pos & offset_mask) + ph_sizeof_dp(sv1, m) >
                    m->pgsize)
                    return PH_SMPGSIZE;
                ph_save_datapoint(new_dp, m);

                Np = 0;
                memcpy(&m->buf[m->file_pos & offset_mask], &Np,
                       sizeof(uint8_t));
                m->file_pos++;
            }
            hfree(sv1->id);
            hfree(sv1->hash);
            ph_free_datapoint(sv1);
        } else { /* put new dp in sv1 pos */
            m->file_pos = start_pos;
            ntype = 0;
            memcpy(&m->buf[m->file_pos & offset_mask], &ntype, sizeof(uint8_t));
            m->file_pos++;
            ph_save_datapoint(new_dp, m);
            ph_save_datapoint(NULL, m);
        }
    } else if (ntype == 1) {
        DP *sv1 = ph_read_datapoint(m);
        DP *sv2 = ph_read_datapoint(m);

        float d1 = hashdist(sv1, new_dp);
        float d2 = hashdist(sv2, new_dp);

        float *M1 = (float *)malloc(LengthM1 * sizeof(float));
        float *M2 = (float *)malloc(LengthM2 * sizeof(float));

        memcpy(M1, &m->buf[m->file_pos & offset_mask],
               LengthM1 * sizeof(float));
        m->file_pos += LengthM1 * sizeof(float);

        memcpy(M2, &m->buf[m->file_pos & offset_mask],
               LengthM2 * sizeof(float));
        m->file_pos += LengthM2 * sizeof(float);

        if (level < m->pathlength) new_dp->path[level] = d1;
        if (level + 1 < m->pathlength) new_dp->path[level + 1] = d2;

        hfree(sv1->id);
        hfree(sv1->hash);
        ph_free_datapoint(sv1);
        hfree(sv2->id);
        hfree(sv2->hash);
        ph_free_datapoint(sv2);

        int pivot1, pivot2;
        uint8_t filenumber;
        off_t child_pos, curr_pos;
        off_t start_pos = m->file_pos;
        MVPRetCode retcode;
        /* check <= each M1 pivot */
        for (pivot1 = 0; pivot1 < LengthM1; pivot1++) {
            if (d1 <= M1[pivot1]) {
                /* check <= each M2 pivot */
                for (pivot2 = 0; pivot2 < LengthM1; pivot2++) {
                    if (d2 <= M2[pivot2 + pivot1 * LengthM1]) {
                        /* determine pos from which to read filenumber and
                         * offset */
                        curr_pos =
                            start_pos + (pivot2 + pivot1 * m->branchfactor) *
                                            (sizeof(uint8_t) + sizeof(off_t));
                        m->file_pos = curr_pos;
                        memcpy(&filenumber, &m->buf[m->file_pos & offset_mask],
                               sizeof(uint8_t));
                        m->file_pos++;
                        memcpy(&child_pos, &m->buf[m->file_pos & offset_mask],
                               sizeof(off_t));
                        m->file_pos += sizeof(off_t);
                        if ((filenumber == 0) && (child_pos == 0)) {
                            FileIndex child_offset;
                            retcode = ph_save_mvptree(m, &new_dp, 1, 0,
                                                      level + 2, &child_offset);
                            if (retcode == PH_SUCCESS) {
                                memcpy(&m->buf[curr_pos++ & offset_mask],
                                       &(child_offset.fileno), sizeof(uint8_t));
                                memcpy(&m->buf[curr_pos & offset_mask],
                                       &(child_offset.offset), sizeof(off_t));
                            }
                            goto addcleanup;
                        } else {
                            /* save position and remap to new file/position */
                            MVPFile m2;
                            retcode = _ph_map_mvpfile(filenumber, child_pos, m,
                                                      &m2, 0);
                            if (ret == PH_SUCCESS) {
                                retcode =
                                    ph_add_mvptree(&m2, new_dp, level + 2);
                                if (retcode != PH_SUCCESS) {
                                    goto addcleanup;
                                }
                                retcode = _ph_unmap_mvpfile(
                                    filenumber, m->file_pos, m, &m2);
                                if (retcode != PH_SUCCESS) {
                                    goto addcleanup;
                                }
                            } else {
                                goto addcleanup;
                            }
                        }
                    }
                }
                /* check > last M2 pivot */
                if (d2 > M2[LengthM1 - 1 + pivot1 * LengthM1]) {
                    curr_pos =
                        start_pos +
                        (m->branchfactor - 1 + pivot1 * m->branchfactor) *
                            (sizeof(uint8_t) + sizeof(off_t));
                    m->file_pos = curr_pos;
                    memcpy(&filenumber, &m->buf[m->file_pos & offset_mask],
                           sizeof(uint8_t));
                    m->file_pos++;
                    memcpy(&child_pos, &m->buf[m->file_pos & offset_mask],
                           sizeof(off_t));
                    m->file_pos += sizeof(off_t);
                    if ((filenumber == 0) && (child_pos == 0)) {
                        FileIndex child_offset;
                        retcode = ph_save_mvptree(m, &new_dp, 1, 0, level + 2,
                                                  &child_offset);
                        if (retcode == PH_SUCCESS) {
                            memcpy(&m->buf[curr_pos++ & offset_mask],
                                   &(child_offset.fileno), sizeof(uint8_t));
                            memcpy(&m->buf[curr_pos & offset_mask],
                                   &(child_offset.offset), sizeof(off_t));
                        }
                        goto addcleanup;
                    } else {
                        MVPFile m2;
                        retcode =
                            _ph_map_mvpfile(filenumber, child_pos, m, &m2, 0);
                        if (retcode == PH_SUCCESS) {
                            retcode = ph_add_mvptree(&m2, new_dp, level + 2);
                            if (retcode != PH_SUCCESS) {
                                goto addcleanup;
                            }
                            retcode = _ph_unmap_mvpfile(filenumber, m->file_pos,
                                                        m, &m2);
                            if (retcode != PH_SUCCESS) {
                                goto addcleanup;
                            }
                        } else {
                            goto addcleanup;
                        }
                    }
                }
            }
        }
        /* check > last M1 pivot */
        if (d1 > M1[LengthM1 - 1]) {
            /*check <= each M2 pivot */
            for (pivot2 = 0; pivot2 < LengthM1; pivot2++) {
                if (d2 <= M2[pivot2 + LengthM1 * LengthM1]) {
                    curr_pos =
                        start_pos + (pivot2 + LengthM1 * m->branchfactor) *
                                        (sizeof(uint8_t) + sizeof(off_t));
                    m->file_pos = curr_pos;
                    memcpy(&filenumber, &m->buf[m->file_pos & offset_mask],
                           sizeof(uint8_t));
                    m->file_pos++;
                    memcpy(&child_pos, &m->buf[m->file_pos & offset_mask],
                           sizeof(off_t));
                    m->file_pos += sizeof(off_t);
                    if ((filenumber == 0) && (child_pos == 0)) {
                        FileIndex child_offset;
                        retcode = ph_save_mvptree(m, &new_dp, 1, 0, level + 2,
                                                  &child_offset);
                        if (retcode == PH_SUCCESS) {
                            memcpy(&m->buf[curr_pos++ & offset_mask],
                                   &(child_offset.fileno), sizeof(uint8_t));
                            memcpy(&m->buf[curr_pos & offset_mask],
                                   &(child_offset.offset), sizeof(off_t));
                        }
                        goto addcleanup;
                    } else {
                        MVPFile m2;
                        retcode =
                            _ph_map_mvpfile(filenumber, child_pos, m, &m2, 0);
                        if (retcode == PH_SUCCESS) {
                            retcode = ph_add_mvptree(&m2, new_dp, level + 2);
                            if (retcode != PH_SUCCESS) {
                                goto addcleanup;
                            }
                            retcode = _ph_unmap_mvpfile(filenumber, m->file_pos,
                                                        m, &m2);
                            if (retcode != PH_SUCCESS) {
                                goto addcleanup;
                            }
                        } else {
                            goto addcleanup;
                        }
                    }
                }
            }

            /* check > last M2 pivot */
            if (d2 > M2[LengthM1 - 1 + LengthM1 * LengthM1]) {
                curr_pos = start_pos +
                           (m->branchfactor - 1 + LengthM1 * m->branchfactor) *
                               (sizeof(uint8_t) + sizeof(off_t));
                m->file_pos = curr_pos;
                memcpy(&filenumber, &m->buf[m->file_pos & offset_mask],
                       sizeof(uint8_t));
                m->file_pos++;
                memcpy(&child_pos, &m->buf[m->file_pos & offset_mask],
                       sizeof(off_t));
                m->file_pos += sizeof(off_t);
                if ((filenumber == 0) && (child_pos == 0)) {
                    FileIndex child_offset;
                    retcode = ph_save_mvptree(m, &new_dp, 1, 0, level + 2,
                                              &child_offset);
                    if (retcode == PH_SUCCESS) {
                        memcpy(&m->buf[curr_pos++ & offset_mask],
                               &(child_offset.fileno), sizeof(uint8_t));
                        memcpy(&m->buf[curr_pos & offset_mask],
                               &(child_offset.offset), sizeof(off_t));
                    }
                    goto addcleanup;
                } else {
                    MVPFile m2;
                    retcode = _ph_map_mvpfile(filenumber, child_pos, m, &m2, 0);
                    if (retcode == PH_SUCCESS) {
                        retcode = ph_add_mvptree(&m2, new_dp, level + 2);
                        if (retcode != PH_SUCCESS) {
                            goto addcleanup;
                        }
                        retcode =
                            _ph_unmap_mvpfile(filenumber, m->file_pos, m, &m2);
                        if (retcode != PH_SUCCESS) {
                            goto addcleanup;
                        }
                    } else {
                        goto addcleanup;
                    }
                }
            }
        }
    addcleanup:
        hfree(M1);
        hfree(M2);
    } else {
        return PH_ERRNTYPE;
    }
    return ret;
}

__declspec(dllexport) MVPRetCode
    ph_add_mvptree(MVPFile *m, DP **points, int nbpoints, int &nbsaved) {
    nbsaved = 0;
    if (!m || !points || nbpoints <= 0) {
        return PH_ERRARG;
    }
    /* open main file */
    char mainfile[MAX_PATH];
    snprintf(mainfile, sizeof(mainfile), "%s.mvp", m->filename);

    m->fh = CreateFile(mainfile, GENERIC_READ | GENERIC_WRITE,
                       FILE_SHARE_READ | FILE_SHARE_WRITE, NULL, OPEN_EXISTING,
                       FILE_ATTRIBUTE_NORMAL, NULL);
    if (m->fh == INVALID_HANDLE_VALUE) {
        return PH_ERRFILEOPEN;
    }
    char *filename = strrchr(m->filename, '\\');
    if (filename == NULL) {
        filename = m->filename;
    }
    char fm_objname[32];
    snprintf(fm_objname, sizeof(fm_objname), "%s%d", filename, 0);
    HANDLE fmhandle =
        CreateFileMapping(m->fh, NULL, PAGE_READWRITE, 0, 0, fm_objname);
    if (fmhandle == NULL) {
        return PH_ERRMMAPCREATE;
    }
    /*use host pg size until pg size of file can be read  */
    m->pgsize = (off_t)getpagesize();

    /* using first file */
    m->filenumber = 0;
    /* map to first page */
    m->file_pos = 0;

    m->buf =
        (char *)MapViewOfFile(fmhandle, FILE_MAP_ALL_ACCESS, 0, 0, m->pgsize);
    if (m->buf == NULL) {
        return PH_ERRMMAPVIEW;
    }
    /* read header within first HeaderSize bytes */
    char tag[17];
    int version;
    ;
    int int_pgsize;
    int leaf_pgsize;
    uint8_t bf, pl, lc, nbdbfiles, type;

    memcpy(tag, &m->buf[m->file_pos], 16);
    tag[16] = '\0';
    m->file_pos += 16;

    memcpy(&version, &m->buf[m->file_pos], sizeof(int));
    m->file_pos += sizeof(int);

    memcpy(&int_pgsize, &m->buf[m->file_pos], sizeof(int));
    m->file_pos += sizeof(int);

    memcpy(&leaf_pgsize, &m->buf[m->file_pos], sizeof(int));
    m->file_pos += sizeof(int);

    off_t nb_pos = m->file_pos;
    memcpy(&nbdbfiles, &m->buf[m->file_pos++], 1);

    memcpy(&bf, &m->buf[m->file_pos++], 1);

    memcpy(&pl, &m->buf[m->file_pos++], 1);

    memcpy(&lc, &m->buf[m->file_pos++], 1);

    memcpy(&type, &m->buf[m->file_pos++], 1);

    m->branchfactor = bf;
    m->pathlength = pl;
    m->leafcapacity = lc;
    m->hash_type = (HashType)type;
    m->nbdbfiles = nbdbfiles;

    m->file_pos = HeaderSize;

    /* add path to datapoints */
    for (int i = 0; i < nbpoints; i++) {
        points[i]->path = (float *)malloc((int)(pl) * sizeof(float));
    }

    /* remap to true pg size used in making file */
    if (!UnmapViewOfFile(m->buf)) {
        return PH_ERRMMAPUNVIEW;
    }
    m->pgsize = (off_t)int_pgsize;
    m->buf =
        (char *)MapViewOfFile(fmhandle, FILE_MAP_ALL_ACCESS, 0, 0, m->pgsize);
    if (m->buf == NULL) return PH_ERRMMAPVIEW;

    MVPRetCode retcode = PH_SUCCESS;

    for (int i = 0; i < nbpoints; i++) {
        retcode = PH_SUCCESS;
        m->file_pos = HeaderSize;
        retcode = ph_add_mvptree(m, points[i], 0);
        if (retcode != PH_SUCCESS) {
            fprintf(stderr, "unable to add %d mvp errorcode %d\n", i, retcode);
            break;
        }
        nbsaved++;
    }

    /* save new nbdbfiles if new files added */
    memcpy(&m->buf[nb_pos], &m->nbdbfiles, 1);

    if (!FlushViewOfFile(m->buf, m->pgsize)) {
        return PH_ERRMMAPFLUSH;
    }

    if (!UnmapViewOfFile(m->buf)) {
        return PH_ERRMMAPUNVIEW;
    }
    if (!CloseHandle(fmhandle)) {
        return PH_ERRFILECLOSE;
    }
    if (!CloseHandle(m->fh)) {
        return PH_ERRFILECLOSE;
    }

    return retcode;
}

__declspec(dllexport) TxtHashPoint *ph_texthash(const char *filename,
                                                int *nbpoints) {
    int count;
    TxtHashPoint *TxtHash = NULL;
    TxtHashPoint WinHash[WindowLength];
    char kgram[KgramLength];

    FILE *pfile = fopen(filename, "r");
    if (!pfile) {
        return NULL;
    }
    struct stat fileinfo;
    fstat(fileno(pfile), &fileinfo);
    count = fileinfo.st_size - WindowLength + 1;
    count = (int)(0.01 * count);
    int d;
    ulong64 hashword = 0ULL;

    TxtHash = (TxtHashPoint *)malloc(count * sizeof(struct ph_hash_point));
    if (!TxtHash) {
        return NULL;
    }
    *nbpoints = 0;
    int i, first = 0;
    int text_index = 0;
    int win_index = 0;
    for (i = 0; i < KgramLength; i++) { /* calc first kgram */
        d = fgetc(pfile);
        if (d == EOF) {
            free(TxtHash);
            return NULL;
        }
        if (d <= 47) /*skip cntrl chars*/
            continue;
        if (((d >= 58) && (d <= 64)) || ((d >= 91) && (d <= 96)) ||
            (d >= 123)) /*skip punct*/
            continue;
        if ((d >= 65) && (d <= 90)) /*convert upper to lower case */
            d = d + 32;

        kgram[i] = (char)d;
        hashword = ROTATELEFT(hashword, delta);
        hashword = hashword ^ textkeys[d];
    }

    WinHash[win_index].hash = hashword;
    WinHash[win_index++].index = text_index;
    struct ph_hash_point minhash;
    minhash.hash = ULLONG_MAX;
    minhash.index = 0;
    struct ph_hash_point prev_minhash;
    prev_minhash.hash = ULLONG_MAX;
    prev_minhash.index = 0;

    while ((d = fgetc(pfile)) != EOF) { /*remaining kgrams */
        text_index++;
        if (d == EOF) {
            free(TxtHash);
            return NULL;
        }
        if (d <= 47) /*skip cntrl chars*/
            continue;
        if (((d >= 58) && (d <= 64)) || ((d >= 91) && (d <= 96)) ||
            (d >= 123)) /*skip punct*/
            continue;
        if ((d >= 65) && (d <= 90)) /*convert upper to lower case */
            d = d + 32;
        hashword = hashword << delta;
        hashword = hashword ^ ((ulong64)d);
        ulong64 oldsym = (ulong64)kgram[first % KgramLength];
        oldsym = oldsym << delta * KgramLength;
        hashword = hashword ^ oldsym;
        kgram[first % KgramLength] = (char)d;
        first++;

        WinHash[win_index % WindowLength].hash = hashword;
        WinHash[win_index % WindowLength].index = text_index;
        win_index++;

        if (win_index >= WindowLength) {
            minhash.hash = ULLONG_MAX;
            for (i = win_index; i < win_index + WindowLength; i++) {
                if (WinHash[i % WindowLength].hash <= minhash.hash) {
                    minhash.hash = WinHash[i % WindowLength].hash;
                    minhash.index = WinHash[i % WindowLength].index;
                }
            }
            if (minhash.hash != prev_minhash.hash) {
                TxtHash[(*nbpoints)].hash = minhash.hash;
                TxtHash[(*nbpoints)++].index = minhash.index;
                prev_minhash.hash = minhash.hash;
                prev_minhash.index = minhash.index;

            } else {
                TxtHash[*nbpoints].hash = prev_minhash.hash;
                TxtHash[(*nbpoints)++].index = prev_minhash.index;
            }
            win_index = 0;
        }
    }

    fclose(pfile);
    return TxtHash;
}

__declspec(dllexport) TxtMatch *ph_compare_text_hashes(TxtHashPoint *hash1,
                                                       int N1,
                                                       TxtHashPoint *hash2,
                                                       int N2, int *nbmatches) {
    int max_matches = (N1 >= N2) ? N1 : N2;
    TxtMatch *found_matches =
        (TxtMatch *)malloc(max_matches * sizeof(TxtMatch));
    if (!found_matches) {
        return NULL;
    }

    *nbmatches = 0;
    int i, j;
    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            if (hash1[i].hash == hash2[j].hash) {
                int m = i + 1;
                int n = j + 1;
                int cnt = 1;
                while ((m < N1) && (n < N2) &&
                       (hash1[m++].hash == hash2[n++].hash)) {
                    cnt++;
                }
                found_matches[*nbmatches].first_index = i;
                found_matches[*nbmatches].second_index = j;
                found_matches[*nbmatches].length = cnt;
                (*nbmatches)++;
            }
        }
    }
    return found_matches;
}

int ph_num_threads() {
    SYSTEM_INFO si;
    GetSystemInfo(&si);
    return si.dwNumberOfProcessors;
}
void ph_mvp_init(MVPFile *m) {
    m->branchfactor = 2;
    m->pathlength = 5;
    m->leafcapacity = 40;
    m->pgsize = 8192; /* use host page size */
    return;
}