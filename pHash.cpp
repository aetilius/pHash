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
#include "cimgffmpeg.h"

const char phash_project[] = "pHash %d.%d.%d. Copyright 2008-2009 David Starkweather & Evan Klinger";
char phash_version[255] = {0};
const char* ph_about(){
	int major, minor, point;
	if(PHASH_VERSION < 100)
		major = 0;
	else
		major = PHASH_VERSION/100;

	if(PHASH_VERSION < 10)
	{
		minor = 0;
		point = PHASH_VERSION;
	}
	else
	{
		minor = (PHASH_VERSION%100)/10;
		point = PHASH_VERSION%10;
	}
	snprintf(phash_version, sizeof(phash_version), phash_project, major, minor, point);
    return phash_version;
}

int ph_radon_projections(const CImg<uint8_t> &img,int N,Projections &projs){

    int width = img.dimx();
    int height = img.dimy();
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
        double theta = k*cimg::valuePI/N;
        double alpha = std::tan(theta);
        for (int x=0;x < D;x++){
	    double y = alpha*(x-x_off);
            int yd = (int)std::floor(y + ROUNDING_FACTOR(y));
            if ((yd + y_off >= 0)&&(yd + y_off < height) && (x < width)){
		ptr_radon_map->at(k,x) = img(x,yd + y_off);
                nb_per_line[k] += 1;
	    }
            if ((yd + x_off >= 0) && (yd + x_off < width) && (k != N/4) && (x < height)){
		ptr_radon_map->at(N/2-k,x) = img(yd + x_off,x);
                nb_per_line[N/2-k] += 1;
	    }
	}
    }
    int j= 0;
    for (int k=3*N/4;k<N;k++){
	double theta = k*cimg::valuePI/N;
        double alpha = std::tan(theta);
        for (int x=0;x < D;x++){
	    double y = alpha*(x-x_off);
            int yd = (int)std::floor(y + ROUNDING_FACTOR(y));
            if ((yd + y_off >= 0)&&(yd + y_off < height) && (x < width)){
		ptr_radon_map->at(k,x) = img(x,yd + y_off);
                nb_per_line[k] += 1;
	    }
            if ((y_off - yd >= 0)&&(y_off - yd<width)&&(2*y_off-x>=0)&&(2*y_off-x<height)&&(k!=3*N/4)){
		ptr_radon_map->at(k-j,x) = img(-yd+y_off,-(x-y_off)+y_off);
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
    int D = projection_map.dimy();

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
	    double temp = R[n]*cos((cimg::valuePI*(2*n+1)*k)/(2*N));
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
	    num  += (x_coeffs[i]-meanx)*(y_coeffs[(i-d)%N]-meany);
            denx += pow((x_coeffs[i]-meanx),2);
            deny += pow((y_coeffs[(i-d)%N]-meany),2);
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

int ph_image_digest(const CImg<uint8_t> &img,double sigma, double gamma,Digest &digest, int N){
    
    int result = EXIT_FAILURE;
    CImg<uint8_t> graysc;
    if (img.dimv() >= 3){
	graysc = img.get_RGBtoYCbCr().channel(0);
    }
    else if (img.dimv() == 1){
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
    
    CImg<uint8_t> *src = new CImg<uint8_t>(file);
    int result = ph_image_digest(*src,sigma,gamma,digest,N);
    delete src;
    return result;

}

int ph_compare_images(const CImg<uint8_t> &imA,const CImg<uint8_t> &imB,double &pcc, double sigma, double gamma,int N,double threshold){

    int result = 0;
    Digest digestA;
    if (ph_image_digest(imA,sigma,gamma,digestA,N) < 0)
	goto cleanup;

    Digest digestB;
    if (ph_image_digest(imB,sigma,gamma,digestB,N) < 0)
	goto cleanup;

    if (ph_crosscorr(digestA,digestB,pcc,threshold) < 0)
	goto cleanup;

    if  (pcc  > threshold)
        result = 1;

cleanup:

    delete &imA;
    delete &imB;
    free(digestA.coeffs);
    free(digestB.coeffs);
    return result;
}

int ph_compare_images(const char *file1, const char *file2,double &pcc, double sigma, double gamma, int N,double threshold){

    CImg<uint8_t> *imA = new CImg<uint8_t>(file1);
    CImg<uint8_t> *imB = new CImg<uint8_t>(file2);
    
    int res = ph_compare_images(*imA,*imB,pcc,sigma,gamma,N,threshold);

    return res;
}

CImg<float>* ph_dct_matrix(const int N){
    CImg<float> *ptr_matrix = new CImg<float>(N,N,1,1,1/sqrt((float)N));
    for (int x=0;x<N;x++){
	for (int y=1;y<N;y++){
	    ptr_matrix->at(x,y) = SQRT_TWO*cos((cimg::valuePI/2/N)*y*(2*(x+1)));
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
    if (src.dimv() >= 3){
        img = src.RGBtoYCbCr().channel(0).get_convolve(meanfilter);
    } else if (img.dimv() ==1){
	img = src.get_convolve(meanfilter);
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

int ph_dct_videohash(const char* file,ulong64 &hash){


    long nb_frames = GetNumberVideoFrames(file);

    if (nb_frames <= 0)
	return -1;

    CImgList<uint8_t> *pframeList = new CImgList<uint8_t>();

    int frames_read = ReadFrames(file,pframeList,0,nb_frames,1,nb_frames);
    if (frames_read < 0){
	printf("error reading frames\n");
	exit(1);
    }

    CImg<float> video_cube = pframeList->get_append('z');
    CImg<float> meanfilter(1,1,41,1,(float)(1.0/41.0));

    video_cube.convolve(meanfilter).resize(32,32,64);

    CImgList<float> *pDCTList = new CImgList<float>();

    CImg<float> *C = ph_dct_matrix(32);
    CImg<float> C_transp = C->get_transpose();
    CImg<float> *C2 = ph_dct_matrix(64);

    cimg_forZ(video_cube,Z){
	CImg<float> current = video_cube.get_slice(Z);
	CImg<float> dctcurrent = (*C)*current*(C_transp);
	(*pDCTList) << dctcurrent;
    }
    video_cube = pDCTList->get_append('z'); 
    pDCTList->clear();
    video_cube.permute_axes("xzyv");
    cimg_forZ(video_cube,Z){
	CImg<float> current = video_cube.get_slice(Z);
	CImg<float> dctcurrent = (*C2)*current;
	(*pDCTList) << dctcurrent;
    }
    video_cube.permute_axes("xzyv");

    CImg<float> coeffs = video_cube.get_crop(1,1,1,4,4,4).unroll('x');
    float median = coeffs.median();
    hash        = 0x0000000000000000ULL;
    ulong64 one = 0x0000000000000001ULL;
    for (int i=0;i<64;i++){
	if (coeffs(i) >= median)
	    hash |= one;
	one = one << 1;
    }

   
    delete pframeList;
    delete pDCTList;

    return 0;
}

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

int ph_rash_videodigest(const char* file,CImg<uint8_t> *p_videodigest){

    int S = 10;
    int L = 50;
    int alpha1 = 3;
    int alpha2 = 2;
    
    CImgList<uint8_t> *pframes = new CImgList<uint8_t>();
    
    long N = GetNumberVideoFrames(file);
    printf("number of video frames %li\n",N);

    int frames_read = ReadFrames(file,pframes,0,N,1,N);
    if (frames_read < 0){
	printf("unable to read frames\n");
	exit(1);
    }
    printf("frames read %d\n",frames_read);
    printf("frame list is size = %d\n",pframes->size);

#ifdef max
#undef max

    CImg<double> dist(N,1,1,1,0);
    CImg<double> prev(64,1,1,1,0);
    cimglist_for(*pframes,I){
	CImg<int> hist = pframes->at(I).get_histogram(64);
        CImg<double> hist_normd = hist/hist.max();
	CImg<int> diff = hist_normd - prev;
	dist(I) = abs(diff.norm(1));
	prev = hist_normd;
    }

    uint8_t bnds[N];
    bnds[0] = 1;
    bnds[N-1] = 1;
    int lstart,lend,gstart,gend;
    for (int k=1;k<N-1;k++){
        lstart = (k-S < 0)?0:k-S;
	lend   = (k+S > N-1)?N-1:k+S;
	gstart = (k-L < 0)?0:k-L;
	gend   = (k+L > N-1)?N-1:k+L;
	CImg<double> local_win = dist.get_crop(lstart,0,lend,0);
	CImg<double> global_win = dist.get_crop(gstart,0,gend,0);
	double Tg = global_win.mean() + alpha1*global_win.variance();
        int local_win_size = local_win.dimx();
	double Tl = alpha2*local_win.kth_smallest(local_win_size-1);
	double localmax = local_win.max();
        double thresh = (Tg >= Tl) ? Tg : Tl;
	if  ((dist(k) >=  thresh) && (dist(k) >= localmax))
	    bnds[k] = 1;
	else
	    bnds[k] = 0;
    }
#endif
    printf("boundaries\n");
    for (int i=0;i<N;i++){
	printf("%d",bnds[i]);
    }
    printf("END\n");
    delete pframes;

    return 0;
}
