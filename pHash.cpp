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

#ifdef HAVE_VIDEO_HASH
#include "cimgffmpeg.h"
#endif

const char phash_project[] = "%s. Copyright 2008-2009 Aetilius, Inc.";
char phash_version[255] = {0};
const char* ph_about(){
	if(phash_version[0] != 0)
		return phash_version;
	
	snprintf(phash_version, sizeof(phash_version), phash_project, PACKAGE_STRING);
	return phash_version;
}
#ifdef HAVE_IMAGE_HASH
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
#endif
#if defined(HAVE_VIDEO_HASH) && defined(HAVE_IMAGE_HASH)
int ph_dct_videohash(const char* file,ulong64 &hash){


    long nb_frames = GetNumberVideoFrames(file);

    if (nb_frames <= 0)
	return -1;

    CImgList<uint8_t> *pframeList = new CImgList<uint8_t>();

    int frames_read = ReadFrames(file,pframeList,0,nb_frames,1,nb_frames);
    if (frames_read < 0){
	delete pframeList;
	return -1;
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

int ph_rash_videodigest(const char* file,CImg<uint8_t> *p_videodigest){

    int S = 10;
    int L = 50;
    int alpha1 = 3;
    int alpha2 = 2;
    
    CImgList<uint8_t> *pframes = new CImgList<uint8_t>();
    
    long N = GetNumberVideoFrames(file);

    int frames_read = ReadFrames(file,pframes,0,N,1,N);
    if (frames_read < 0){
	delete pframes;
	return -1;
    }

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
    delete pframes;

    return 0;
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

DP* ph_malloc_datapoint(){
    DP* dp = (DP*)malloc(sizeof(DP));
    dp->hash = 0LLU;
    dp->id = NULL;
    dp->path = (int*)malloc(PathLength*sizeof(int));
    for (int i=0;i<PathLength;i++){
	dp->path[i] = -1;
    }
    return dp;
}
void ph_free_datapoint(DP *dp){
    if (dp){
	free(dp->path);
	free(dp->id);
	free(dp);
    }
    return;
}


Node* ph_malloc_leaf(){
    Node *newnode = (Node*)malloc(sizeof(Node));
    newnode->leaf.node_type = 0;
    return (newnode);
}

Node* ph_malloc_inner(){
    Node *newnode = (Node*)malloc(sizeof(Node));
    newnode->internal.node_type = 1;
    newnode->internal.M1 = (float*)malloc(LengthM1*sizeof(float));
    newnode->internal.M2 = (float*)malloc(LengthM2*sizeof(float));
    return (newnode);
}

void ph_free_node(Node *node){
    if (node->leaf.node_type == 0){
        ph_free_datapoint(node->leaf.sv1);
	ph_free_datapoint(node->leaf.sv2);
	free(node->leaf.D1);
	free(node->leaf.D2);
	free(node->leaf.points);
    } else if (node->internal.node_type == 1){
        ph_free_datapoint(node->leaf.sv1);
	ph_free_datapoint(node->leaf.sv2);
	for (int i=0;i<node->internal.Nc;i++){
	    Node *child_node = (Node*)node->internal.child[i];
	    ph_free_node(child_node);
	}
	free(node->internal.M1);
	free(node->internal.M2);
	free(node->internal.child);
    }
    free(node);
}


Node* ph_create_MVPtree(DP **points, int nbpoints, int level=0){
    Node *new_node = NULL;
    int Np = nbpoints - 2;
    if (nbpoints <= 0){
	return NULL;
    }
    if (nbpoints <= K+2) { //create leaf node
        new_node = ph_malloc_leaf();
	if (Np >= 0){
	     new_node->leaf.sv1 = points[0];
	     new_node->leaf.points = (DP**)malloc(Np*sizeof(DP**));
	     new_node->leaf.D1 = (int*)calloc(Np,sizeof(int));
	     new_node->leaf.D2 = (int*)calloc(Np,sizeof(int));
	     new_node->leaf.Np = Np;
	     int max_dist = 0;
	     int max_pos = 1;
	     int d;
	     for (int i=0;i<nbpoints;i++){
		d = ph_hamming_distance(new_node->leaf.sv1->hash, points[i]->hash);
		if (d > max_dist){
		    max_dist = d;
		    max_pos = i;
		}
	    }
	    new_node->leaf.sv2 = *(points+max_pos);
	    int index = 1;
	    for (int i=0;i<Np;i++){
		if (index == max_pos)
		    index++;
		d = ph_hamming_distance(new_node->leaf.sv1->hash, points[index]->hash);
		new_node->leaf.D1[i] = d;
		d = ph_hamming_distance(new_node->leaf.sv2->hash, points[index]->hash);
		new_node->leaf.D2[i] = d;
		new_node->leaf.points[i] = points[index];
		index++;
	    }
	} else if (Np == -1){ //Np == -1
	    new_node->leaf.sv1 = points[0];
	    new_node->leaf.sv2 = NULL;
	    new_node->leaf.points = NULL;
	    new_node->leaf.D1 = NULL;
	    new_node->leaf.D2 = NULL;
	    new_node->leaf.Np = 0;
	} 
    }else { //create internal_node - nodes greater than K 
	new_node = ph_malloc_inner();
	//chose sv1,sv2
	new_node->internal.sv1 = points[0];
	int max_dist = 0, min_dist = 100,d;
	int max_pos = 1;
        int *dist = (int*)malloc(nbpoints*sizeof(int));
	for (int i=1;i<nbpoints;i++){
	    dist[i] = ph_hamming_distance(new_node->internal.sv1->hash,points[i]->hash);
            if (level <= PathLength)
		points[i]->path[level] = dist[i];
	    if ((dist[i] < min_dist)&&(dist[i] != 0))
		min_dist = dist[i];
	    if (dist[i] > max_dist){
		max_dist = dist[i];
		max_pos = i;
	    }
	}
	new_node->internal.sv2 = *(points + max_pos);
	for (int i=0;i<nbpoints;i++){
	    if (level < PathLength){
		d = ph_hamming_distance(new_node->internal.sv2->hash,points[i]->hash);
		points[i]->path[level+1] = d;
	    }
	}
	//determine median values
	float slice = (max_dist - min_dist)/BranchFactor;
	float step = slice;
	
	for (int i=0;i < LengthM1;i++){
	    new_node->internal.M1[i] = min_dist + slice;
	    slice += step;
	}

	//set up 1st bins
	DP ***bins = (DP***)malloc(BranchFactor*sizeof(DP***));
	int *mlens = (int*)calloc(BranchFactor,sizeof(int));
	for (int i=0;i < BranchFactor;i++){
	    bins[i] = (DP**)malloc(((int)Np)*sizeof(DP**));
	}

        //sort into bins
	for (int i=1;i<nbpoints;i++){
	    if (i == max_pos)
		continue;
	    int cur_dist = dist[i];
	    for (int j=0;j < LengthM1;j++){
		if (cur_dist <= new_node->internal.M1[j]){
		    bins[j][mlens[j]] = points[i];
		    mlens[j]++;
                    break;
		}
	    }
	    if (cur_dist > new_node->internal.M1[LengthM1-1]){
		bins[BranchFactor-1][mlens[BranchFactor-1]] = points[i];
		mlens[BranchFactor-1]++;
	    }
	}

	//set up 2nd set of bins
	DP ***bins2 = (DP***)malloc(BranchFactor*sizeof(DP***));
	int *mlens2 = (int*)calloc(BranchFactor,sizeof(int));
	for (int i=0;i < BranchFactor;i++){
	    bins2[i] = (DP**)malloc(Np*sizeof(DP**));
	}

        new_node->internal.child = (Node**)malloc(Fanout*sizeof(Node*));
        new_node->internal.Nc = Fanout;

	//for each row of bin
        for (int i=0;i < BranchFactor;i++){
	    int d_len = mlens[i];
	    for (int j=0;j < BranchFactor;j++){ //reset bins2 array 
		mlens2[j] = 0;
	    }
	    dist = (int*)realloc(dist,d_len*sizeof(int));
	    max_dist = 0;
	    min_dist = 100;
	    for (int j=0;j<d_len;j++){
	       dist[j]=ph_hamming_distance(new_node->internal.sv2->hash,bins[i][j]->hash);
	       if (dist[j] > max_dist)
		   max_dist = dist[j];
	       if (dist[j] < min_dist)
		   min_dist = dist[j];
	    }
	    slice = (max_dist - min_dist)/BranchFactor;
	    step = slice;
	    for (int j=0;j < LengthM1;j++){
		new_node->internal.M2[j + i*(LengthM1)] = slice + min_dist;
		slice += step;
	    }

	    //sort bins[i][]  into bins2 
	    for (int j=0;j < d_len;j++){
	        DP *current = bins[i][j];
		for (int k=0;k < LengthM1 ;k++){
		    if (dist[j] <= new_node->internal.M2[k + i*(LengthM1)]){
			bins2[k][mlens2[k]] = current;
			mlens2[k]++;
			break;
		    }
		}
	   	if (dist[j] > new_node->internal.M2[(LengthM1-1)+i*(LengthM1)]){
		    bins2[BranchFactor-1][mlens2[BranchFactor-1]] = current;
		    mlens2[BranchFactor-1]++;
		}
	    } 
	    //create child nodes
	    for (int j=0;j < BranchFactor;j++){
		Node *child_node = ph_createMVPtree(bins2[j],mlens2[j],level+2);
	        new_node->internal.child[j+i*BranchFactor] = child_node;
	    }
	}
	//cleanup
	free(dist);
	free(mlens);
	free(mlens2);
	free(bins);
	free(bins2);
    }
    return new_node;
}


Node* ph_addDPtoMVPtree(Node *tree,DP *dp,int level){
    Node *ret_node = NULL;
    if (tree->leaf.node_type ==0){
	int Np = tree->leaf.Np;
	DP *sv1 = tree->leaf.sv1;
	DP *sv2 = tree->leaf.sv2;
	DP **dplist = (DP**)malloc((Np+3)*sizeof(DP**));

	int count = 0;
	dplist[count++] = sv1;
        dplist[count++] = sv2;
	for (int i=0;i<Np;i++){
	    dplist[count++] = tree->leaf.points[i];
	}
	dplist[count++] = dp;
        ret_node = ph_create_MVPtree(dplist,count,level+2);
        free(dplist);
    } else {
	int d1 = ph_hamming_distance(tree->internal.sv1->hash,dp->hash);
	int d2 = ph_hamming_distance(tree->internal.sv2->hash,dp->hash);
	if (level <= PathLength){
	    dp->path[level] = d1;
	}
	if (level < PathLength){
	    dp->path[level+1] = d2;
	}
	for (int pivot1=0;pivot1<LengthM1;pivot1++){
	    if (d1 <= tree->internal.M1[pivot1]){
		for (int pivot2=0;pivot2<LengthM1;pivot2++){
		    if (d2 <= tree->internal.M2[pivot2 + pivot1*(LengthM1)]){
			Node *subtree = tree->internal.child[pivot2 + pivot1*BranchFactor];
			Node *newsub = ph_addDPtoMVPtree(subtree,dp,level+2);
			free(subtree);
			tree->internal.child[pivot2+pivot1*BranchFactor] = newsub;
		    }
		}
                if (d2 > tree->internal.M2[LengthM1-1 + pivot1*(LengthM1)]){
		    Node *subtree = tree->internal.child[BranchFactor-1+pivot1*BranchFactor];
		    Node *newsub = ph_addDPtoMVPtree(subtree,dp, level+2);
		    free(subtree);
		    tree->internal.child[BranchFactor-1+pivot1*BranchFactor] = newsub;
		}
	    }
	}
	if (d1 > tree->internal.M1[LengthM1-1]){
	    for (int pivot2=0;pivot2<LengthM1;pivot2++){
		if (d2 <= tree->internal.M2[pivot2 + (LengthM1)*(LengthM1)]){
		    Node *subtree = tree->internal.child[pivot2 + (LengthM1)*BranchFactor];
		    Node *newsub = ph_addDPtoMVPtree(subtree,dp,level+2);
		    free(subtree);
		    tree->internal.child[pivot2+(LengthM1)*BranchFactor] = newsub;
		}
	    }
	    if (d2 > tree->internal.M2[LengthM1-1+(LengthM1)*(LengthM1)]){
		Node *subtree = tree->internal.child[BranchFactor-1+(LengthM1)*BranchFactor];
		Node *newsub = ph_addDPtoMVPtree(subtree,dp,level+2);
		free(subtree);
		tree->internal.child[BranchFactor-1+(LengthM1)*BranchFactor] = newsub;
	    }
	}
	ret_node = tree;
    }
    return ret_node;
}

DP** ph_read_imagehashes(const char *dirname,int capacity, int &count){

    count = 0;
    struct dirent *dir_entry;
    DIR *dir = opendir(dirname);
    if (!dir)
	exit(1);

    DP **hashlist = (DP**)malloc(capacity*sizeof(DP**));
    if (!hashlist)
	exit(1);

    DP *dp = NULL;
    errno = 0;
    ulong64 tmphash = 0;
    char path[100];
    path[0] = '\0';

    while ((dir_entry = readdir(dir)) != 0){
	if (strcmp(dir_entry->d_name,".") && strcmp(dir_entry->d_name,"..")){
	    strcat(path, dirname);
	    strcat(path, "/");
	    strcat(path, dir_entry->d_name);
	    if (ph_dct_imagehash(path, tmphash) < 0)  //calculate the hash
		continue;
	    dp = ph_malloc_datapoint();
	    dp->id = strdup(path);
	    dp->hash = tmphash;
	    hashlist[count] = dp;
	    count++;
	}
	errno = 0;
        path[0]='\0';
    }
    if (errno)
	exit(1);

    return hashlist;

}

char** ph_readfilenames(const char *dirname,int cap,int &count){
    count = 0;
    struct dirent *dir_entry;
    DIR *dir = opendir(dirname);
    if (!dir)
	exit(1);

    char **filelist = (char**)malloc(cap*sizeof(char*));
    if (!filelist)
	exit(1);

    errno = 0;
    char path[100];
    path[0] = '\0';

    while ((dir_entry = readdir(dir)) != 0){
	if (strcmp(dir_entry->d_name,".") && strcmp(dir_entry->d_name,"..")){
	    strcat(path, dirname);
	    strcat(path, "/");
	    strcat(path, dir_entry->d_name);
	    filelist[count++] = strdup(path);
	}
	errno = 0;
        path[0]='\0';
    }
    if (errno)
	exit(1);

    return filelist;
}


int ph_saveDP(DP *dp, FILE *pfile){
    int length_p = P;
    char c = '\0';
    if (!pfile)
	return -1;
    if (!dp){
	fwrite(&c,sizeof(char),1,pfile);
	return 0;
    }
    fputs(dp->id,pfile);
    fputc(c,pfile);
    fwrite(&(dp->hash), sizeof(ulong64), 1, pfile);
    fwrite(&length_p, sizeof(int), 1, pfile);
    fwrite(dp->path, sizeof(int),length_p, pfile);
    return 0;
}
DP * ph_readDP(FILE *pfile){
    if (!pfile)
	return NULL;
    char id[100];
    ulong64 hash;
    int length_p;
    DP *dp = ph_malloc_datapoint();
    char c;
    int i=0;
    do {
	fread(&c, sizeof(char),1,pfile);
	id[i++]=c;
    } while (c != '\0');
    if (id[0] == '\0')
	return NULL;
    dp->id = strdup(id);
    fread(&hash, sizeof(ulong64),1,pfile);
    dp->hash = hash;
    fread(&length_p, sizeof(int),1,pfile);
    fread(dp->path, sizeof(int),length_p,pfile);
   
    return dp;
}


Node* ph_readMVPtree(FILE *pfile){
    Node *new_node = NULL;
    if (!pfile){
	return NULL;
    }
    int node_type;
    fread(&node_type,sizeof(int),1,pfile);
    if (node_type == 0){ //create leaf
	new_node = ph_malloc_leaf();
	int Np;
	new_node->leaf.sv1 = ph_readDP(pfile);
	new_node->leaf.sv2 = ph_readDP(pfile);
	fread(&Np, sizeof(int),1,pfile);
	new_node->leaf.Np = Np;
	new_node->leaf.points = (DP**)malloc(Np*sizeof(DP*));
	new_node->leaf.D1 = (int*)malloc(Np*sizeof(int));
	new_node->leaf.D2 = (int*)malloc(Np*sizeof(int));
	for (int i=0;i<Np;i++){
	    new_node->leaf.points[i] = ph_readDP(pfile);
	}
	fread(new_node->leaf.D1,sizeof(int), Np, pfile);
	fread(new_node->leaf.D2,sizeof(int), Np, pfile);
    } else { //creat internal node
	new_node = ph_malloc_inner();
	new_node->internal.sv1 = ph_readDP(pfile);
	new_node->internal.sv2 = ph_readDP(pfile);

	fread(new_node->internal.M1, sizeof(float),LengthM1,pfile);
	fread(new_node->internal.M2, sizeof(float),LengthM2,pfile);
	new_node->internal.Nc = Fanout;
	new_node->internal.child = (Node**)malloc(Fanout*sizeof(Node*));
	for (int i=0;i<Fanout;i++){
	    new_node->internal.child[i] = ph_readMVPtree(pfile);
	}
    }
    return new_node;
}


Node* ph_readMVPtree(const char *filename){
    FILE *pfile = fopen(filename,"r");
    if (!pfile)
	return NULL;
    fread(&BranchFactor,sizeof(int),1,pfile);
    fread(&PathLength,sizeof(int),1,pfile);
    fread(&LeafCapacity,sizeof(int),1,pfile);
    LengthM1 = BranchFactor - 1;
    LengthM2 = BranchFactor*LengthM1;
    Fanout   = BranchFactor*BranchFactor;
    Node *node = ph_readMVPtree(pfile);
    fclose(pfile);
    return (node);
}

int ph_saveMVPtree(Node *tree, FILE *pfile){
    if (!pfile){
	return -1;
    }
    if (!tree){
	return 0;
    }
    
    if (tree->leaf.node_type==0){
	fwrite(&(tree->leaf.node_type),sizeof(int),1,pfile);
	int Np = tree->leaf.Np;
	if (ph_saveDP(tree->leaf.sv1,pfile) < 0)
	    return -1;
	if (ph_saveDP(tree->leaf.sv2,pfile) < 0)
	    return -1;
	fwrite(&Np, sizeof(int),1,pfile);
	for (int i=0;i<Np;i++){
	    DP *current = tree->leaf.points[i];
	    if (ph_saveDP(current,pfile) < 0)
		return -1;
	}
        fwrite(tree->leaf.D1, sizeof(int), Np, pfile);
	fwrite(tree->leaf.D2, sizeof(int), Np, pfile);
    } else {
	fwrite(&(tree->internal.node_type),sizeof(int),1,pfile);
	if (ph_saveDP(tree->internal.sv1, pfile) < 0)
	    return -1;
	if (ph_saveDP(tree->internal.sv2, pfile) < 0)
	    return -1;
	fwrite(tree->internal.M1, sizeof(float),LengthM1, pfile);
	fwrite(tree->internal.M2, sizeof(float),LengthM2, pfile);
	for (int i=0;i < Fanout;i++){
	    Node *current = tree->internal.child[i];
	    if (ph_saveMVPtree(current, pfile) < 0)
		return -1;
	}

    }
    return 0;
}


int ph_saveMVPtree(Node *tree, const char *filename){
    int m = M;
    int p = P;
    int k = K;
    FILE *pfile = fopen(filename,"w");
    if (!pfile)
	return -1;
    fwrite(&m,sizeof(int),1,pfile);
    fwrite(&p,sizeof(int),1,pfile);
    fwrite(&k,sizeof(int),1,pfile);
    int result = ph_saveMVPtree(tree,pfile);
    fclose(pfile);
    return result;
}

void ph_printMVPtree(Node *tree){
    if (!tree){
	printf("\nNULL\n");
	return;
    }
    if (tree->leaf.node_type == 0){ //if leaf node
	printf("\nleaf node\n");
	DP *sv1=NULL, *sv2=NULL;
	sv1 = tree->leaf.sv1;
	sv2 = tree->leaf.sv2;
	if (sv1){
	    printf("sv1: %s\n",sv1->id);
	}
	if (sv2){
	    printf("sv2: %s\n",sv2->id);
	}
	printf("number points: Np = %d\n", tree->leaf.Np);
        for (int i=0;i< tree->leaf.Np;i++){
	    printf("point:%s,D[%d]=%d, D[%d]=%d\n", tree->leaf.points[i]->id,i,tree->leaf.D1[i],i,tree->leaf.D2[i]);
	}
    } else {
	printf("\ninternal node\n");
	DP *sv1 = tree->internal.sv1;
	DP *sv2 = tree->internal.sv2;
	if (sv1){
	    printf("sv1: %s\n", sv1->id);
	}
	if (sv2){
	    printf("sv2: %s\n", sv2->id);
	}
	printf("First level pivot points:\n");
	for (int i=0;i<LengthM1;i++){
	    printf(" M1[%d]=%f ",i,tree->internal.M1[i]);
	}
	printf("\n");
	printf("Second level pivot points:\n");
	for (int i=0;i<LengthM2;i++){
	    printf(" M2[%d]=%f ",i,tree->internal.M2[i]);
	}
	printf("\n");

	for (int i=0;i < Fanout;i++){
	    printf(" %d child node", i);
	    Node *child = (Node*)tree->internal.child[i];
	    ph_printMVPtree(child);
	}
    }
}


int ph_queryMVPtree(Node *tree,DP *query,int r,int k,int *path,DP **results,int &count,int level){
    int d1,d2;
    if (level == 0){
	count = 0;
    }
    if (!results){
	return -1;
    }
    if (!tree){
	return 0;
    }
    if (tree->leaf.node_type == 0){ //if leaf node
	DP *sv1=NULL, *sv2=NULL;
	sv1 = tree->leaf.sv1;
	sv2 = tree->leaf.sv2;
	d1 = ph_hamming_distance(query->hash,sv1->hash);
	if (d1 <= r){
	    results[count++]=sv1;
	}
	if (sv2){
	    d2 = ph_hamming_distance(query->hash,sv2->hash);
	    if (d2 <= r){
		    results[count++]=sv2;
	    }
	    for (int i=0;i < tree->leaf.Np;i++){
		DP *current = tree->leaf.points[i];
		int da = tree->leaf.D1[i];
		int db = tree->leaf.D2[i];
		int include = 1;
		if (((d1-r <= da)&&(d1+r >= da)) && ((d2 -r <= db)&&(d2+r >= db))){
		    for (int j=0;j < PathLength;j++){
			if ((path[j]-r <= current->path[j]) && (path[j]+r >= current->path[j]))
			    if (ph_hamming_distance(query->hash, current->hash)>r){
				include = 0;
				break;
			    }
		    }
		    if (include){
		    results[count++] = current;
		    }
		}
	    }
	}
    } else {
	DP *sv1 = tree->internal.sv1;
	DP *sv2 = tree->internal.sv2;
	d1 = ph_hamming_distance(query->hash, sv1->hash);
	d2 = ph_hamming_distance(query->hash, sv2->hash);
	if (d1 <= r){
	    results[count++]=sv1;
	}
	if (d2 <= r){
	    results[count++]=sv2;
        }
        if (level < PathLength)
	    path[level]=d1;
	if (level < PathLength-1)
	    path[level+1] = d2;

	int pivot1;
	int pivot2;

	//check <= each M1 pivot
	for (pivot1=0;pivot1 < LengthM1;pivot1++){
	    if (d1 - r <= tree->internal.M1[pivot1]){
		//check <= each M2 pivot
		for (pivot2=0;pivot2 < LengthM1;pivot2++){
		    if (d2 - r <= tree->internal.M2[pivot2+pivot1*(LengthM1)]){
			Node *child = (Node*)tree->internal.child[pivot2+pivot1*BranchFactor];
			ph_queryMVPtree(child,query,r,k,path,results,count, level+2);
		    }
		}
		//check > last M2 pivot
		if (d2+r >= tree->internal.M2[LengthM1-1+pivot1*(LengthM1)]){
		    Node *child = (Node*)tree->internal.child[BranchFactor-1+pivot1*BranchFactor];
		    ph_queryMVPtree(child,query,r,k,path,results,count,level+2);
		}
	    }
	}

        if (d1+r >= tree->internal.M1[LengthM1-1])  {

	    //check > last M1 pivot
	    for (pivot2=0;pivot2<LengthM1;pivot2++){
		if (d2-r <= tree->internal.M2[pivot2+(LengthM1)*(LengthM1)]){
		    Node *child = (Node*)tree->internal.child[pivot2+(LengthM1)*BranchFactor];
		    ph_queryMVPtree(child,query,r,k,path,results,count,level+2);
		}
	        
	    }
	    if (d2+r > tree->internal.M2[LengthM1-1+(LengthM1)*(LengthM1)]){
		Node *child = (Node*)tree->internal.child[BranchFactor-1+(LengthM1)*BranchFactor];
		ph_queryMVPtree(child,query,r,k,path,results,count,level+2);
	    }
	}

    }

    return count;
}
