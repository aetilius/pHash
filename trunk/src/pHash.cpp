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
    D Grant Starkweather - dstarkweather@phash.org

*/

#include "pHash.h"
#include "config.h"
#ifdef HAVE_VIDEO_HASH
#include "cimgffmpeg.h"
#endif

#ifdef HAVE_PTHREAD
#include <pthread.h>

int ph_num_threads()
{
	int numCPU = 1;
#ifdef __GLIBC__
		numCPU = sysconf( _SC_NPROCESSORS_ONLN );
#else
		int mib[2];
		size_t len; 

		mib[0] = CTL_HW;
		mib[1] = HW_AVAILCPU;

		sysctl(mib, 2, &numCPU, &len, NULL, 0);

		if( numCPU < 1 ) 
		{
     			mib[1] = HW_NCPU;
     			sysctl( mib, 2, &numCPU, &len, NULL, 0 );

		     	if( numCPU < 1 )
     			{
          			numCPU = 1;
     			}
		}

#endif
	return numCPU;
}
#endif

const char phash_project[] = "%s. Copyright 2008-2010 Aetilius, Inc.";
char phash_version[255] = {0};
const char* ph_about(){
	if(phash_version[0] != 0)
		return phash_version;
	
	snprintf(phash_version, sizeof(phash_version), phash_project, PACKAGE_STRING);
	return phash_version;
}
#ifdef HAVE_IMAGE_HASH
int ph_radon_projections(const CImg<uint8_t> &img,int N,Projections &projs){

    int width = img.width();
    int height = img.height();
    int D = (width > height)?width:height;
    float x_center = (float)width/2;
    float y_center = (float)height/2;
    int x_off = (int)std::floor(x_center + ROUNDING_FACTOR(x_center));
    int y_off = (int)std::floor(y_center + ROUNDING_FACTOR(y_center));

    projs.R = new CImg<uint8_t>(N,D,1,1,0);
    projs.nb_pix_perline = (int*)calloc(N,sizeof(int));

    if (!projs.R || !projs.nb_pix_perline)
	return EXIT_FAILURE;

    projs.size = N;

    CImg<uint8_t> *ptr_radon_map = projs.R;
    int *nb_per_line = projs.nb_pix_perline;

    for (int k=0;k<N/4+1;k++){
        double theta = k*cimg::PI/N;
        double alpha = std::tan(theta);
        for (int x=0;x < D;x++){
	    double y = alpha*(x-x_off);
            int yd = (int)std::floor(y + ROUNDING_FACTOR(y));
            if ((yd + y_off >= 0)&&(yd + y_off < height) && (x < width)){
		*ptr_radon_map->data(k,x) = img(x,yd + y_off);
                nb_per_line[k] += 1;
	    }
            if ((yd + x_off >= 0) && (yd + x_off < width) && (k != N/4) && (x < height)){
		*ptr_radon_map->data(N/2-k,x) = img(yd + x_off,x);
                nb_per_line[N/2-k] += 1;
	    }
	}
    }
    int j= 0;
    for (int k=3*N/4;k<N;k++){
	double theta = k*cimg::PI/N;
        double alpha = std::tan(theta);
        for (int x=0;x < D;x++){
	    double y = alpha*(x-x_off);
            int yd = (int)std::floor(y + ROUNDING_FACTOR(y));
            if ((yd + y_off >= 0)&&(yd + y_off < height) && (x < width)){
		*ptr_radon_map->data(k,x) = img(x,yd + y_off);
                nb_per_line[k] += 1;
	    }
            if ((y_off - yd >= 0)&&(y_off - yd<width)&&(2*y_off-x>=0)&&(2*y_off-x<height)&&(k!=3*N/4)){
		*ptr_radon_map->data(k-j,x) = img(-yd+y_off,-(x-y_off)+y_off);
                nb_per_line[k-j] += 1;
	    }
            
	}
        j += 2;
    }

    return EXIT_SUCCESS;

}
int ph_feature_vector(const Projections &projs, Features &fv)
{

    CImg<uint8_t> *ptr_map = projs.R;
    CImg<uint8_t> projection_map = *ptr_map;
    int *nb_perline = projs.nb_pix_perline;
    int N = projs.size;
    int D = projection_map.height();

    fv.features = (double*)malloc(N*sizeof(double));
    fv.size = N;
    if (!fv.features)
	return EXIT_FAILURE;

    double *feat_v = fv.features;
    double sum = 0.0;
    double sum_sqd = 0.0;
    for (int k=0; k < N; k++){
	double line_sum = 0.0;
        double line_sum_sqd = 0.0;
        int nb_pixels = nb_perline[k];
	for (int i=0;i<D;i++){
	    line_sum += projection_map(k,i);
    	    line_sum_sqd += projection_map(k,i)*projection_map(k,i);
	}
	feat_v[k] = (line_sum_sqd/nb_pixels) - (line_sum*line_sum)/(nb_pixels*nb_pixels);
        sum += feat_v[k];
        sum_sqd += feat_v[k]*feat_v[k];
    }
    double mean = sum/N;
    double var  = sqrt((sum_sqd/N) - (sum*sum)/(N*N));

    for (int i=0;i<N;i++){
    	feat_v[i] = (feat_v[i] - mean)/var;
    }

    return EXIT_SUCCESS;
} 
int ph_dct(const Features &fv,Digest &digest)
{
    int N = fv.size;
    const int nb_coeffs = 40;

    digest.coeffs = (uint8_t*)malloc(nb_coeffs*sizeof(uint8_t));
    if (!digest.coeffs)
	return EXIT_FAILURE;

    digest.size = nb_coeffs;

    double *R = fv.features;

    uint8_t *D = digest.coeffs;

    double D_temp[nb_coeffs];
    double max = 0.0;
    double min = 0.0;
    for (int k = 0;k<nb_coeffs;k++){
	double sum = 0.0;
        for (int n=0;n<N;n++){
	    double temp = R[n]*cos((cimg::PI*(2*n+1)*k)/(2*N));
            sum += temp;
	}
        if (k == 0)
	    D_temp[k] = sum/sqrt((double)N);
        else
            D_temp[k] = sum*SQRT_TWO/sqrt((double)N);
        if (D_temp[k] > max)
            max = D_temp[k];
        if (D_temp[k] < min)
            min = D_temp[k];
    }
       
    for (int i=0;i<nb_coeffs;i++){

	D[i] = (uint8_t)(UCHAR_MAX*(D_temp[i] - min)/(max - min));

    }
    
    return EXIT_SUCCESS;
}

int ph_crosscorr(const Digest &x,const Digest &y,double &pcc,double threshold){

    int N = y.size;
    int result = 0;

    uint8_t *x_coeffs = x.coeffs;
    uint8_t *y_coeffs = y.coeffs;

    double *r = new double[N];
    double sumx = 0.0;
    double sumy = 0.0;
    for (int i=0;i < N;i++){
	sumx += x_coeffs[i];
        sumy += y_coeffs[i];
    }
    double meanx = sumx/N;
    double meany = sumy/N;
    double max = 0;
    for (int d=0;d<N;d++){
        double num = 0.0;
        double denx = 0.0;
        double deny = 0.0;
	for (int i=0;i<N;i++){
	    num  += (x_coeffs[i]-meanx)*(y_coeffs[(N+i-d)%N]-meany);
            denx += pow((x_coeffs[i]-meanx),2);
            deny += pow((y_coeffs[(N+i-d)%N]-meany),2);
	}
        r[d] = num/sqrt(denx*deny);
        if (r[d] > max)
	    max = r[d];
    }
    delete[] r;
    pcc = max;
    if (max > threshold)
	    result = 1;

    return result;
}

#ifdef max
#undef max
#endif

int _ph_image_digest(const CImg<uint8_t> &img,double sigma, double gamma,Digest &digest, int N){
    
    int result = EXIT_FAILURE;
    CImg<uint8_t> graysc;
    if (img.spectrum() >= 3){
	graysc = img.get_RGBtoYCbCr().channel(0);
    }
    else if (img.spectrum() == 1){
	graysc = img;
    }
    else {
	return result;
    }
	
 
    graysc.blur((float)sigma);
 
    (graysc/graysc.max()).pow(gamma);
     
    Projections projs;
    if (ph_radon_projections(graysc,N,projs) < 0)
	goto cleanup;
 
    Features features;
    if (ph_feature_vector(projs,features) < 0)
	goto cleanup;
    
    if (ph_dct(features,digest) < 0)
        goto cleanup;
 
    result = EXIT_SUCCESS;

cleanup:
    free(projs.nb_pix_perline);
    free(features.features);

    delete projs.R;
    return result;
}

#define max(a,b) (((a)>(b))?(a):(b))

int ph_image_digest(const char *file, double sigma, double gamma, Digest &digest, int N){
    
    CImg<uint8_t> src(file);
	int res = -1;
    		int result = _ph_image_digest(src,sigma,gamma,digest,N);
    		res = result;
	return res;
}

int _ph_compare_images(const CImg<uint8_t> &imA,const CImg<uint8_t> &imB,double &pcc, double sigma, double gamma,int N,double threshold){

    int result = 0;
    Digest digestA;
    if (_ph_image_digest(imA,sigma,gamma,digestA,N) < 0)
	goto cleanup;

    Digest digestB;
    if (_ph_image_digest(imB,sigma,gamma,digestB,N) < 0)
	goto cleanup;

    if (ph_crosscorr(digestA,digestB,pcc,threshold) < 0)
	goto cleanup;

    if  (pcc  > threshold)
        result = 1;

cleanup:

    free(digestA.coeffs);
    free(digestB.coeffs);
    return result;
}

int ph_compare_images(const char *file1, const char *file2,double &pcc, double sigma, double gamma, int N,double threshold){

    CImg<uint8_t> *imA = new CImg<uint8_t>(file1);
    CImg<uint8_t> *imB = new CImg<uint8_t>(file2);
    
    int res = _ph_compare_images(*imA,*imB,pcc,sigma,gamma,N,threshold);

    delete imA;
    delete imB;
    return res;
}

CImg<float>* ph_dct_matrix(const int N){
    CImg<float> *ptr_matrix = new CImg<float>(N,N,1,1,1/sqrt((float)N));
    const float c1 = sqrt(2.0/N); 
    for (int x=0;x<N;x++){
	for (int y=1;y<N;y++){
	    *ptr_matrix->data(x,y) = c1*cos((cimg::PI/2/N)*y*(2*x+1));
	}
    }
    return ptr_matrix;
}

int ph_dct_imagehash(const char* file,ulong64 &hash){

    if (!file){
	return -1;
    }
    CImg<uint8_t> src;
    try {
	src.load(file);
    } catch (CImgIOException ex){
	return -1;
    }
    CImg<float> meanfilter(7,7,1,1,1);
    CImg<float> img;
    if (src.spectrum() == 3){
        img = src.RGBtoYCbCr().channel(0).get_convolve(meanfilter);
    } else if (src.spectrum() == 4){
	int width = img.width();
        int height = img.height();
        int depth = img.depth();
	img = src.crop(0,0,0,0,width-1,height-1,depth-1,2).RGBtoYCbCr().channel(0).get_convolve(meanfilter);
    } else {
	img = src.channel(0).get_convolve(meanfilter);
    }

    img.resize(32,32);
    CImg<float> *C  = ph_dct_matrix(32);
    CImg<float> Ctransp = C->get_transpose();

    CImg<float> dctImage = (*C)*img*Ctransp;

    CImg<float> subsec = dctImage.crop(1,1,8,8).unroll('x');;
   
    float median = subsec.median();
    ulong64 one = 0x0000000000000001;
    hash = 0x0000000000000000;
    for (int i=0;i< 64;i++){
	float current = subsec(i);
        if (current > median)
	    hash |= one;
	one = one << 1;
    }
  
    delete C;

    return 0;
}

#ifdef HAVE_PTHREAD
void *ph_image_thread(void *p)
{
        slice *s = (slice *)p;
        for(int i = 0; i < s->n; ++i)
        {
                DP *dp = (DP *)s->hash_p[i];
		ulong64 hash;
		int ret = ph_dct_imagehash(dp->id, hash);
                dp->hash = (ulong64*)malloc(sizeof(hash));
		memcpy(dp->hash, &hash, sizeof(hash));
                dp->hash_length = 1;
        }
}

DP** ph_dct_image_hashes(char *files[], int count, int threads)
{
       	if(!files || count <= 0)
                return NULL;

        int num_threads;
        if(threads > count)
        {
                num_threads = count;
        }
	else if(threads > 0)
        {
                num_threads = threads;
        }
	else
	{
                num_threads = ph_num_threads();
        }

	DP **hashes = (DP**)malloc(count*sizeof(DP*));

        for(int i = 0; i < count; ++i)
        {
                hashes[i] = (DP *)malloc(sizeof(DP));
                hashes[i]->id = strdup(files[i]);
	}

	pthread_t thds[num_threads];

        int rem = count % num_threads;
        int start = 0;
        int off = 0;
        slice *s = new slice[num_threads];
        for(int n = 0; n < num_threads; ++n)
        {
                off = (int)floor((count/(float)num_threads) + (rem>0?num_threads-(count % num_threads):0));

                s[n].hash_p = &hashes[start];
                s[n].n = off;
                s[n].hash_params = NULL;
                start += off;
                --rem;
                pthread_create(&thds[n], NULL, ph_image_thread, &s[n]);
        }
	for(int i = 0; i < num_threads; ++i)
        {
                pthread_join(thds[i], NULL);
        }
	delete[] s;

        return hashes;

}
#endif


#endif

#if defined(HAVE_VIDEO_HASH) && defined(HAVE_IMAGE_HASH)


CImgList<uint8_t>* ph_getKeyFramesFromVideo(const char *filename){

    long N =  GetNumberVideoFrames(filename);

    if (N < 0){
	return NULL;
    }
    
    float frames_per_sec = 0.5*fps(filename);
    if (frames_per_sec < 0){
	return NULL;
    }

    int step = (int)(frames_per_sec + ROUNDING_FACTOR(frames_per_sec));
    long nbframes = (long)(N/step);

    float *dist = (float*)malloc((nbframes)*sizeof(float));
    if (!dist){
	return NULL;
    }
    CImg<float> prev(64,1,1,1,0);

    VFInfo st_info;
    st_info.filename = filename;
    st_info.nb_retrieval = 100;
    st_info.step = step;
    st_info.pixelformat = 0;
    st_info.pFormatCtx = NULL;
    st_info.width = -1;
    st_info.height = -1;

    CImgList<uint8_t> *pframelist = new CImgList<uint8_t>();
    if (!pframelist){
	return NULL;
    }
    int nbread = 0;
    int k=0;
    do {
	nbread = NextFrames(&st_info, pframelist);
	if (nbread < 0){
	    delete pframelist;
	    free(dist);
	    return NULL;
	}
	unsigned int i = 0;
        while ((i < pframelist->size()) && (k < nbframes)){
	    CImg<uint8_t> current = pframelist->at(i++);
	    CImg<float> hist = current.get_histogram(64,0,255);
            float d = 0.0;
	    dist[k] = 0.0;
	    cimg_forX(hist,X){
		d =  hist(X) - prev(X);
		d = (d>=0) ? d : -d;
		dist[k] += d;
		prev(X) = hist(X);
	    }
            k++;
	}
        pframelist->clear();
    } while ((nbread >= st_info.nb_retrieval)&&(k < nbframes));
    vfinfo_close(&st_info);

    int S = 10;
    int L = 50;
    int alpha1 = 3;
    int alpha2 = 2;
    int s_begin, s_end;
    int l_begin, l_end;
    uint8_t *bnds = (uint8_t*)malloc(nbframes*sizeof(uint8_t));
    if (!bnds){
	delete pframelist;
	free(dist);
	return NULL;
    }

    int nbboundaries = 0;
    k = 1;
    bnds[0] = 1;
    do {
	s_begin = (k-S >= 0) ? k-S : 0;
	s_end   = (k+S < nbframes) ? k+S : nbframes-1;
	l_begin = (k-L >= 0) ? k-L : 0;
	l_end   = (k+L < nbframes) ? k+L : nbframes-1;

	/* get global average */
	float ave_global, sum_global = 0.0, dev_global = 0.0;
	for (int i=l_begin;i<=l_end;i++){
	    sum_global += dist[i];
	}
	ave_global = sum_global/((float)(l_end-l_begin+1));

	/*get global deviation */
	for (int i=l_begin;i<=l_end;i++){
	    float dev = ave_global - dist[i];
	    dev = (dev >= 0) ? dev : -1*dev;
	    dev_global += dev;
	}
	dev_global = dev_global/((float)(l_end-l_begin+1));

	/* global threshold */
	float T_global = ave_global + alpha1*dev_global;

	/* get local maximum */
	int localmaxpos = s_begin;
	for (int i=s_begin;i<=s_end;i++){
	    if (dist[i] > dist[localmaxpos])
		localmaxpos = i;
	}
        /* get 2nd local maximum */
	int localmaxpos2 = s_begin;
	float localmax2 = 0;
	for (int i=s_begin;i<=s_end;i++){
	    if (i == localmaxpos)
		continue;
	    if (dist[i] > localmax2){
		localmaxpos2 = i;
		localmax2 = dist[i];
	    }
	}
        float T_local = alpha2*dist[localmaxpos2];
	float Thresh = (T_global >= T_local) ? T_global : T_local;

	if ((dist[k] == dist[localmaxpos])&&(dist[k] > Thresh)){
	    bnds[k] = 1;
	    nbboundaries++;
	}
	else {
	    bnds[k] = 0;
	}
	k++;
    } while ( k < nbframes-1);
    bnds[nbframes-1]=1;
    nbboundaries += 2;

    int start = 0;
    int end = 0;
    int nbselectedframes = 0;
    do {
	/* find next boundary */
	do {end++;} while ((bnds[end]!=1)&&(end < nbframes));
       
	/* find min disparity within bounds */
	int minpos = start+1;
	for (int i=start+1; i < end;i++){
	    if (dist[i] < dist[minpos])
		minpos = i;
	}
	bnds[minpos] = 2;
	nbselectedframes++;
	start = end;
    } while (start < nbframes-1);

    st_info.nb_retrieval = 1;
    st_info.width = 32;
    st_info.height = 32;
    k = 0;
    do {
	if (bnds[k]==2){
	    if (ReadFrames(&st_info, pframelist, k*st_info.step,k*st_info.step + 1) < 0){
		delete pframelist;
		free(dist);
		return NULL;
	    }
	}
	k++;
    } while (k < nbframes);
    vfinfo_close(&st_info);

    free(bnds);
    bnds = NULL;
    free(dist);
    dist = NULL;

    return pframelist;
}


ulong64* ph_dct_videohash(const char *filename, int &Length){

    CImgList<uint8_t> *keyframes = ph_getKeyFramesFromVideo(filename);
    if (keyframes == NULL)
	return NULL;

    Length = keyframes->size();

    ulong64 *hash = (ulong64*)malloc(sizeof(ulong64)*Length);
    CImg<float> *C = ph_dct_matrix(32);
    CImg<float> Ctransp = C->get_transpose();
    CImg<float> dctImage;
    CImg<float> subsec;
    CImg<uint8_t> currentframe;

    for (unsigned int i=0;i < keyframes->size(); i++){
	currentframe = keyframes->at(i);
	currentframe.blur(1.0);
	dctImage = (*C)*(currentframe)*Ctransp;
	subsec = dctImage.crop(1,1,8,8).unroll('x');
	float med = subsec.median();
	hash[i] =     0x0000000000000000;
	ulong64 one = 0x0000000000000001;
	for (int j=0;j<64;j++){
	    if (subsec(j) > med)
		hash[i] |= one;
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

#ifdef HAVE_PTHREAD
void *ph_video_thread(void *p)
{
        slice *s = (slice *)p;
        for(int i = 0; i < s->n; ++i)
        {
                DP *dp = (DP *)s->hash_p[i];
		int N;
		ulong64 *hash = ph_dct_videohash(dp->id, N);
		if(hash)
		{
                	dp->hash = hash;
	                dp->hash_length = N;
		}
		else
		{
			dp->hash = NULL;
			dp->hash_length = 0;
		}
        }
}

DP** ph_dct_video_hashes(char *files[], int count, int threads)
{
       	if(!files || count <= 0)
                return NULL;

        int num_threads;
        if(threads > count)
        {
                num_threads = count;
        }
	else if(threads > 0)
        {
                num_threads = threads;
        }
	else
	{
                num_threads = ph_num_threads();
        }

	DP **hashes = (DP**)malloc(count*sizeof(DP*));

        for(int i = 0; i < count; ++i)
        {
                hashes[i] = (DP *)malloc(sizeof(DP));
                hashes[i]->id = strdup(files[i]);
	}

	pthread_t thds[num_threads];

        int rem = count % num_threads;
        int start = 0;
        int off = 0;
        slice *s = new slice[num_threads];
        for(int n = 0; n < num_threads; ++n)
        {
                off = (int)floor((count/(float)num_threads) + (rem>0?num_threads-(count % num_threads):0));

                s[n].hash_p = &hashes[start];
                s[n].n = off;
                s[n].hash_params = NULL;
                start += off;
                --rem;
                pthread_create(&thds[n], NULL, ph_video_thread, &s[n]);
        }
	for(int i = 0; i < num_threads; ++i)
        {
                pthread_join(thds[i], NULL);
        }
	delete[] s;

        return hashes;

}
#endif

double ph_dct_videohash_dist(ulong64 *hashA, int N1, ulong64 *hashB, int N2, int threshold){

    int den = (N1 <= N2) ? N1 : N2;
    int C[N1+1][N2+1];

    for (int i=0;i<N1+1;i++){
	C[i][0] = 0;
    }
    for (int j=0;j<N2+1;j++){
	C[0][j] = 0;
    }
    for (int i=1;i<N1+1;i++){
	for (int j=1;j<N2+1;j++){
	    int d = ph_hamming_distance(hashA[i-1],hashB[j-1]);
	    if (d <= threshold){
		C[i][j] = C[i-1][j-1] + 1;
	    } else {
		C[i][j] = ((C[i-1][j] >= C[i][j-1])) ? C[i-1][j] : C[i][j-1];
	    }
	}
    }

    double result = (double)(C[N1][N2])/(double)(den);

    return result;
}

#endif

int ph_hamming_distance(const ulong64 hash1,const ulong64 hash2){
    ulong64 x = hash1^hash2;
    const ulong64 m1  = 0x5555555555555555ULL;
    const ulong64 m2  = 0x3333333333333333ULL;
    const ulong64 h01 = 0x0101010101010101ULL;
    const ulong64 m4  = 0x0f0f0f0f0f0f0f0fULL;
    x -= (x >> 1) & m1;
    x = (x & m2) + ((x >> 2) & m2);
    x = (x + (x >> 4)) & m4;
    return (x * h01)>>56;
}

DP* ph_malloc_datapoint(int hashtype){
    DP* dp = (DP*)malloc(sizeof(DP));
    dp->hash = NULL;
    dp->id = NULL;
    dp->hash_type = hashtype;
    return dp;
}
void ph_free_datapoint(DP *dp){
    if (!dp){
	return;
    }
    free(dp);
    dp=NULL;
    return;
}


#ifdef HAVE_IMAGE_HASH

CImg<float>* GetMHKernel(float alpha, float level){
    int sigma = (int)4*pow((float)alpha,(float)level);
    static CImg<float> *pkernel = NULL;
    float xpos, ypos, A;
    if (!pkernel){
	pkernel = new CImg<float>(2*sigma+1,2*sigma+1,1,1,0);
        cimg_forXY(*pkernel,X,Y){
	    xpos = pow(alpha,-level)*(X-sigma);
            ypos = pow(alpha,-level)*(Y-sigma);
            A = xpos*xpos + ypos*ypos;
            pkernel->atXY(X,Y) = (2-A)*exp(-A/2);
	}
    }
    return pkernel;
}

uint8_t* ph_mh_imagehash(const char *filename, int &N,float alpha, float lvl){
    if (filename == NULL){
	return NULL;
    }
    uint8_t *hash = (unsigned char*)malloc(72*sizeof(uint8_t));
    N = 72;

    CImg<uint8_t> src(filename);
    CImg<uint8_t> img;

    if (src.spectrum() == 3){
	img = src.get_RGBtoYCbCr().channel(0).blur(1.0).resize(512,512,1,1,5).get_equalize(256);
    } else{
	img = src.channel(0).get_blur(1.0).resize(512,512,1,1,5).get_equalize(256);
    }
    src.clear();

    CImg<float> *pkernel = GetMHKernel(alpha,lvl);
    CImg<float> fresp =  img.get_correlate(*pkernel);
    img.clear();
    fresp.normalize(0,1.0);
    CImg<float> blocks(31,31,1,1,0);
    for (int rindex=0;rindex < 31;rindex++){
		for (int cindex=0;cindex < 31;cindex++){
			blocks(rindex,cindex) = fresp.get_crop(rindex*16,cindex*16,rindex*16+16-1,cindex*16+16-1).sum();
		}
    }
    int hash_index;
    int nb_ones = 0, nb_zeros = 0;
    int bit_index = 0;
    unsigned char hashbyte = 0;
    for (int rindex=0;rindex < 31-2;rindex+=4){
		CImg<float> subsec;
		for (int cindex=0;cindex < 31-2;cindex+=4){
			subsec = blocks.get_crop(cindex,rindex, cindex+2, rindex+2).unroll('x');
			float ave = subsec.mean();
			cimg_forX(subsec, I){
				hashbyte <<= 1;
				if (subsec(I) > ave){
					hashbyte |= 0x01;
					nb_ones++;
				} else {
					nb_zeros++;
				}
				bit_index++;
				if ((bit_index%8) == 0){
					hash_index = (int)(bit_index/8) - 1; 
					hash[hash_index] = hashbyte;
					hashbyte = 0x00;
				}
			}
		}
	}

    return hash;
}
#endif

char** ph_readfilenames(const char *dirname,int &count){
    count = 0;
    struct dirent *dir_entry;
    DIR *dir = opendir(dirname);
    if (!dir)
        return NULL;

    /*count files */
    while ((dir_entry = readdir(dir)) != NULL){
	if (strcmp(dir_entry->d_name, ".") && strcmp(dir_entry->d_name,".."))
	    count++;
    }
    
    /* alloc list of files */
    char **files = (char**)malloc(count*sizeof(*files));
    if (!files)
	return NULL;

    errno = 0;
    int index = 0;
    char path[1024];
    path[0] = '\0';
    rewinddir(dir);
    while ((dir_entry = readdir(dir)) != 0){
	if (strcmp(dir_entry->d_name,".") && strcmp(dir_entry->d_name,"..")){
	    strcat(path, dirname);
	    strcat(path, "/");
	    strcat(path, dir_entry->d_name);
	    files[index++] = strdup(path);
	}
        path[0]='\0';
    }
    if (errno)
	return NULL;
    closedir(dir);
    return files;
}


int ph_bitcount8(uint8_t val){
    int num = 0;
    while (val){
	++num;
	val &= val - 1;
    }
    return num;
}



double ph_hammingdistance2(uint8_t *hashA, int lenA, uint8_t *hashB, int lenB){
    if (lenA != lenB){
	return -1.0;
    }
    if ((hashA == NULL) || (hashB == NULL) || (lenA <= 0)){
	return -1.0;
    }
    double dist = 0;
    uint8_t D = 0;
    for (int i=0;i<lenA;i++){
	D = hashA[i]^hashB[i];
	dist = dist + (double)ph_bitcount8(D);
    }
    double bits = (double)lenA*8;
    return dist/bits;

}


TxtHashPoint* ph_texthash(const char *filename,int *nbpoints){
    int count;
    TxtHashPoint *TxtHash = NULL;
    TxtHashPoint WinHash[WindowLength];
    char kgram[KgramLength];

    FILE *pfile = fopen(filename,"r");
    if (!pfile){
	return NULL;
    }
    struct stat fileinfo;
    fstat(fileno(pfile),&fileinfo);
    count = fileinfo.st_size - WindowLength + 1;
    count = (int)(0.01*count);
    int d;
    ulong64 hashword = 0ULL;
    
    TxtHash = (TxtHashPoint*)malloc(count*sizeof(struct ph_hash_point));
    if (!TxtHash){
	return NULL;
    }
    *nbpoints=0;
    int i, first=0, last=KgramLength-1;
    int text_index = 0;
    int win_index = 0;
    for (i=0;i < KgramLength;i++){    /* calc first kgram */
	d = fgetc(pfile);
	if (d == EOF){
	    free(TxtHash);
	    return NULL;
	}
	if (d <= 47)         /*skip cntrl chars*/
	    continue;
	if ( ((d >= 58)&&(d <= 64)) || ((d >= 91)&&(d <= 96)) || (d >= 123) ) /*skip punct*/
	    continue;
	if ((d >= 65)&&(d<=90))       /*convert upper to lower case */
	    d = d + 32;
      
	kgram[i] = (char)d;
        hashword = hashword << delta;   /* rotate left or shift left ??? */
        hashword = hashword^textkeys[d];/* right now, rotate breaks it */
    }

    WinHash[win_index].hash = hashword;
    WinHash[win_index++].index = text_index;
    struct ph_hash_point minhash;
    minhash.hash = ULLONG_MAX;
    minhash.index = 0;
    struct ph_hash_point prev_minhash;
    prev_minhash.hash = ULLONG_MAX;
    prev_minhash.index = 0;

    while ((d=fgetc(pfile)) != EOF){    /*remaining kgrams */
        text_index++;
	if (d == EOF){
	    free(TxtHash);
	    return NULL;
	}
	if (d <= 47)         /*skip cntrl chars*/
	    continue;
	if ( ((d >= 58)&&(d <= 64)) || ((d >= 91)&&(d <= 96)) || (d >= 123) ) /*skip punct*/
	    continue;
	if ((d >= 65)&&(d<=90))       /*convert upper to lower case */
	    d = d + 32;

	ulong64 oldsym = textkeys[kgram[first%KgramLength]];

	/* rotate or left shift ??? */
	/* right now, rotate breaks it */
	oldsym = oldsym << delta*KgramLength;
	hashword = hashword << delta;
	hashword = hashword^textkeys[d];
	hashword = hashword^oldsym;
	kgram[last%KgramLength] = (char)d;
	first++;
	last++;

        WinHash[win_index%WindowLength].hash = hashword;
	WinHash[win_index%WindowLength].index = text_index;
	win_index++;

        if (win_index >= WindowLength){
	    minhash.hash = ULLONG_MAX;
	    for (i=win_index;i<win_index+WindowLength;i++){
		if (WinHash[i%WindowLength].hash <= minhash.hash){
		    minhash.hash = WinHash[i%WindowLength].hash;
		    minhash.index = WinHash[i%WindowLength].index;
		}
	    }
            if (minhash.hash != prev_minhash.hash){	 
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

TxtMatch* ph_compare_text_hashes(TxtHashPoint *hash1, int N1, TxtHashPoint *hash2,int N2, int *nbmatches){

    int max_matches = (N1 >= N2) ? N1:N2;
    TxtMatch *found_matches = (TxtMatch*)malloc(max_matches*sizeof(TxtMatch));
    if (!found_matches){
	return NULL;
    }

    *nbmatches = 0;
    int i,j;
    for (i=0;i<N1;i++){
	for (j=0;j<N2;j++){
	    if (hash1[i].hash == hash2[j].hash){
		int m = i + 1;
		int n = j + 1;
                int cnt = 1;
		while((m < N1)&&(n < N2)&&(hash1[m++].hash == hash2[n++].hash)){
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

