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
#include <cmath>
#include "pHash.h"
#ifdef USE_VIDEO_HASH
#include "cimgffmpeg.h"
#endif

const char phash_project[] = "%s. Copyright 2008-2010 Aetilius, Inc.";

char phash_version[255] = {0};

const char* ph_about(){
	if(phash_version[0] != 0)
		return phash_version;
	
	snprintf(phash_version, sizeof(phash_version), phash_project, PACKAGE_STRING);
	return phash_version;
}

#if defined USE_IMAGE_HASH || defined USE_VIDEO_HASH

int ph_hamming_distance(const ulong64 hash1,const ulong64 hash2){
    ulong64 x = hash1^hash2;
    return __builtin_popcountll(x);
}

static CImg<float>* ph_dct_matrix(const int N){
    CImg<float> *ptr_matrix = new CImg<float>(N,N,1,1,1/sqrt((float)N));
    const float c1 = sqrt(2.0/N); 
    for (int x=0;x<N;x++){
		for (int y=1;y<N;y++){
			*ptr_matrix->data(x,y) = c1*cos((cimg::PI/2/N)*y*(2*x+1));
		}
    }
    return ptr_matrix;
}

#endif /* USE_IMAGE_HASH || USE_VIDEO_HASH */

#ifdef USE_IMAGE_HASH

void ph_free_feature(Features &fv){
	delete[] fv.features;
}

void ph_free_digest(Digest &digest){
	delete[] digest.coeffs;
}


int ph_feature_vector(const char *file, const int n_angles, const double sigma, const double gamma, Features &fv){
	const int preset_width = 256;
	const int preset_height = 256;
	const int strip_width = 5;
	
	CImg<float> src;
	try {
		src.load(file);
	} catch (CImgException &ex){
		return -1;
	}

	if (src.spectrum() >= 3){
		src.RGBtoYCbCr().channel(0);
    } else if (src.spectrum() == 1){
		src.channel(0);
    }
    else {
		return -1;
    }

	src.normalize(0,1);
    src.blur((float)sigma, false, true);
	src.pow(gamma);
	src.resize(preset_width, preset_height, 1, 1, 5);

	float x_center = round((float)src.width()/2.0);
	float y_center = round((float)src.width()/2.0);
	
    fv.features = new double[n_angles]; 
    fv.n_features = n_angles;
	if (fv.features == NULL)
		return -1;
	
	/* construct mask of a thin horizontal strip across the center */
	CImg<bool> mask(src.width(), src.height(), 1, 1);
	cimg_forXY(mask, x, y){
		mask(x,y) = false;
		if ((y >= y_center - strip_width) && (y <= y_center + strip_width)){
				mask(x,y) = true;
		} 
	}

	/* rotate mask for each angle */
	for (int i=0;i < n_angles;i++){
		float angle = -(float)i*180.0/(float)n_angles;
		/* rotate mask - 2 for cubic interpolation */
		CImg<uint8_t> rotd_mask = mask.get_rotate(angle, (float)x_center, (float)y_center, 2, 0);

		/* compute variance for pixels under mask*/
		int n_pixels = 0;
		double sum = 0, sumsq = 0;
		cimg_forXY(rotd_mask, x, y){
			if (rotd_mask(x,y)){
				sum = sum + src(x,y);
				sumsq = sumsq + src(x,y)*src(x,y);
				n_pixels += 1;
			}
		}
		fv.features[i] = (sumsq - (sum*sum)/n_pixels)/n_pixels;
	}
	
	return 0;
}

int ph_dct(const Features &fv,Digest &digest){
	const int n_coeffs = 40;
	double sqrt_two = sqrt(2.0);
	double sqrt_n_features = sqrt(fv.n_features);
	
    digest.coeffs = new uint8_t[n_coeffs];
	digest.n_coeffs = n_coeffs;
	if (digest.coeffs == NULL)
		return -1;

	double max_c = numeric_limits<double>::min();
    double min_c = numeric_limits<double>::max();
    
  	double tmp_coeffs[n_coeffs];

	/* calculate dct coefficients, coeffs[k] */
    for (int k = 0;k < n_coeffs;k++){
		double sum = 0.0;
        for (int n=0;n < fv.n_features;n++){
			sum +=  fv.features[n]*cos((cimg::PI*(2*n+1)*k)/(2*fv.n_features));
		}
        if (k == 0)
			tmp_coeffs[k] = sum/sqrt((double)fv.n_features);
        else
            tmp_coeffs[k] = sum*sqrt_two/sqrt_n_features;
        if (tmp_coeffs[k] > max_c)
            max_c = tmp_coeffs[k];
        if (tmp_coeffs[k] < min_c)
            min_c = tmp_coeffs[k];
    }
       
	/* quantize coefficients into bytes */
    for (int i=0;i<n_coeffs;i++){
		digest.coeffs[i] = static_cast<uint8_t>(UCHAR_MAX*(tmp_coeffs[i] - min_c)/(max_c - min_c));
    }
    
    return 0;
}

int ph_image_digest(const char *file,const double sigma, const double gamma,Digest &digest, const int n_angles){

	Features features;
    if (ph_feature_vector(file, n_angles, sigma, gamma, features) < 0)
		return -1;

    int rc = ph_dct(features,digest);
	if (rc < 0){
		ph_free_feature(features);
		return -1;
	}
 
	ph_free_feature(features);
	
    return 0;
}

int ph_peakcrosscorr(const Digest &x,const Digest &y,double &pcc){
	if (x.n_coeffs != y.n_coeffs || x.coeffs == NULL || y.coeffs == NULL)
		return -1;

    uint8_t *x_coeffs = x.coeffs;
    uint8_t *y_coeffs = y.coeffs;
	int N = x.n_coeffs;

	/* calculate mean of two series */
	double mx = 0, my = 0;
	for (int i=0;i < N;i++){
		mx += static_cast<double>(x_coeffs[i])/255.0;
		my += static_cast<double>(y_coeffs[i])/255.0;
		
	}
	mx /= static_cast<double>(N);
	my /= static_cast<double>(N);

	/* calculate dens */
	double sx = 0, sy = 0, den;
	for (int i=0;i < N;i++){
		sx += ((double)x_coeffs[i]/255.0 - mx)*((double)x_coeffs[i]/255.0 - mx);
		sy += ((double)y_coeffs[i]/255.0 - my)*((double)y_coeffs[i]/255.0 - my);
	}
	den = sqrt(sx)*sqrt(sy);

	/* calculate correlation series */
	double cc, max_cc = 0;
	for (int d=0;d < N;d++){
		double sxy = 0;
		for (int i=0;i < N;i++){
			sxy += ((double)x_coeffs[i]/255.0 - mx)*((double)y_coeffs[(i-d+N)%N]/255.0 - my);
		}
		cc = sxy / den;
		if (cc > max_cc)
			max_cc = cc;
	}
	pcc = max_cc;

	return 0;
}

int ph_dct_imagehash(const char* file,ulong64 &hash){

    if (!file)
		return -1;
    
    CImg<uint8_t> src;
    try {
		src.load(file);
    } catch (CImgIOException& ex){
		return -1;
    }
    CImg<float> meanfilter(7,7,1,1,1);
    CImg<float> img;
    if (src.spectrum() == 3){
        img = src.RGBtoYCbCr().channel(0).get_convolve(meanfilter);
    } else if (src.spectrum() == 4){
		int width = src.width();
        int height = src.height();
		img = src.crop(0,0,0,0,width-1,height-1,0,2).RGBtoYCbCr().channel(0).get_convolve(meanfilter);
    } else {
		img = src.channel(0).get_convolve(meanfilter);
    }

    img.resize(32,32);
    CImg<float> *C  = ph_dct_matrix(32);
    CImg<float> Ctransp = C->get_transpose();

    CImg<float> dctImage = (*C)*img*Ctransp;

    CImg<float> subsec = dctImage.crop(1,1,8,8).unroll('x');;
   
    float median = subsec.median();
    hash = 0;
    for (int i=0;i < 64;i++, hash <<= 1){
		float current = subsec(i);
        if (current > median)
			hash |= 0x01;
    }
  
    delete C;

    return 0;
}


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

int ph_bitcount8(uint8_t val){
    int num = 0;
    while (val){
		++num;
		val &= val - 1;
    }
    return num;
}



int ph_hammingdistance2(uint8_t *hashA, int lenA, uint8_t *hashB, int lenB){
    if (lenA != lenB)
		return -1.0;
    
    if ((hashA == NULL) || (hashB == NULL) || (lenA <= 0))
		return -1.0;
    
    int dist = 0;
	for (int i=0;i<lenA;i++){
		uint8_t xord = hashA[i]^hashB[i];
		int d  = ph_bitcount8(xord);
		dist += d;
    }

    return dist;
}
#endif /* USE_IMAGE_HASH */


#ifdef USE_VIDEO_HASH

CImgList<uint8_t>* ph_getKeyFramesFromVideo(const char *filename){

    long N =  GetNumberVideoFrames(filename);

    if (N < 0){
		return NULL;
    }
    
    float frames_per_sec = 0.5*fps(filename);
    if (frames_per_sec < 0){
		return NULL;
    }

    int step =static_cast<int>(round(frames_per_sec)); 
    long nbframes = static_cast<long>(N/step);

	float *dist =  new float[nbframes];
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
			delete [] dist;
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
	uint8_t *bnds = new uint8_t[nbframes];
    if (!bnds){
		delete pframelist;
		delete [] dist;
		return NULL;
    }

    int nbboundaries = 0;
    k = 1;
    bnds[0] = 1;
    do {
		s_begin = (k-S >= 0) ? k-S : 0;
		s_end   = (k+S < nbframes) ? k+S : nbframes-1;
		l_begin=  (k-L >= 0) ? k-L : 0;
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
		int minpos=  start+1;
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
				delete [] dist;
				return NULL;
			}
		}
		k++;
    } while (k < nbframes);
    vfinfo_close(&st_info);

    delete [] bnds;
    bnds = NULL;
    delete [] dist;
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

#endif /* USE_VIDEO_HASH */

#ifdef USE_TEXT_HASH

TxtHashPoint* ph_texthash(const char *filename,int *nbpoints){
    int count;
    TxtHashPoint *TxtHash = NULL;
    TxtHashPoint WinHash[WindowLength];
    char kgram[KgramLength];

    FILE *pfile = fopen(filename,"r");
    if (!pfile)
		return NULL;
    
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
				TxtHash[*nbpoints].index = prev_minhash.index;
			}
			win_index = 0;
		}
    }

    fclose(pfile);
    return TxtHash;
}

TxtMatch* ph_compare_text_hashes(TxtHashPoint *hash1, int N1, TxtHashPoint *hash2,int N2, int *nbmatches){

    int max_matches = (N1 >= N2) ? N1 : N2;
    TxtMatch *found_matches = (TxtMatch*)malloc(max_matches*sizeof(TxtMatch));
    if (!found_matches)
		return NULL;

    *nbmatches = 0;

    for (int i=0;i<N1;i++){
		for (int j=0;j<N2;j++){
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

#endif /* USE_TEXT_HASH */


