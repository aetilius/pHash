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
	X[k]     = S[k] + P[k];
	X[k+N/2] = S[k] - P[k];
    }
    for (k=0;k<N/2;k++){
	X[k] =     S[k] + P[k];
	X[k+N/2] = S[k] - P[k];
    }
    

}


int fft(double *x, int N, complex double *X){

    complex double *twiddle_factors = (complex double*)malloc(sizeof(complex double)*(N/2));
    complex double *Xt = (complex double*)malloc(sizeof(complex double)*N);
    int k;
    for (k=0;k<N/2;k++){
	twiddle_factors[k] = polar_to_complex(1.0, 2.0*PI*k/N);
    }
    fft_calc(N, x, X, Xt, 1, twiddle_factors);

    free(twiddle_factors);
    free(Xt);

    return 0;

}
