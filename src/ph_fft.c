#include "ph_fft.h"

complex double polar_to_complex(double r, double theta){
    complex double result;
    result = r*cos(theta) + r*sin(theta)*I;

    return result;
}

void fft_calc(int N,double *x,complex double *X,complex double *P,int step,complex double *twids){
    complex double *S = P + N/2;
    if (N == 1){
	X[0] = x[0];
	return;
    }
    
    fft_calc(N/2, x,      S,   X,2*step, twids);
    fft_calc(N/2, x+step, P,   X,2*step, twids);	    

    int k;
    for (k=0;k<N/2;k++){
	P[k] = P[k]*twids[k*step];
        nbmults += 2;
    }
    for (k=0;k<N/2;k++){
	X[k] =     S[k] + P[k];
	X[k+N/2] = S[k] - P[k];
	nbadds += 2;
    }
    

}


int fft(double *x, int N, complex double *X){
    if (sizeof(X) != N*sizeof(complex double)){
	return -1;
    }
    complex double *twiddle_factors = (complex double*)malloc(sizeof(complex double)*(N/2));
    complex double *Xt = (complex double*)malloc(sizeof(complex double)*N);
    int k;
    for (k=0;k<N/2;k++){
	twiddle_factors[k] = polar_to_complex(1.0, 2.0*PI*k/N);
    }
    nbadds = 0;
    nbmults = 0;
    fft_calc(N, x, X, Xt, 1, twiddle_factors);

    free(twiddle_factors);
    free(Xt);

    return 0;

}
