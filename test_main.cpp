#include <cstdio>
#include "CImg.h"
#include "pHash.h"

using namespace cimg_library;
int main(int argc, char **argv){
    cimg_usage("pHash Radial Hash projection algorithm");
    const char *msg = ph_about();
    puts(msg);

    const char *filename1 = cimg_option("-f1",(char *)NULL,"first image file name");
    const char *filename2 = cimg_option("-f2",(char *)NULL,"second image file name");
    const int nb_angles  = cimg_option("-n",180,"number of angles");
    double sigma = cimg_option("-s",3.5,"sigma parameter for blurring");
    double gamma = cimg_option("-g",1.0,"gamma parameter for gamma-correction");
    double thresh = cimg_option("-t",0.90,"pcc threshold value");

    printf("filename1: %s\n",filename1);
    printf("filename2: %s\n",filename2);
    printf("nb angles %d\n",nb_angles);
    printf("sigma = %f\n",sigma);
    printf("gamma = %f\n",gamma);
    printf("threshold = %f\n",thresh);

    double pcc_value;
    int result = ph_compare_images(filename1,filename2,pcc_value,sigma,gamma,nb_angles,thresh);

    if (result < 0){
	printf("unable to complete comparison\n");
	exit(1);
    }
    if (result){
	printf("same image. pcc = %f",pcc_value);
        return 1;
    }
    printf("different images.  pcc = %f\n",pcc_value);
    return 0;
}
