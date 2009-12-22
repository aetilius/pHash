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
#include "ph_fft.h"

Complexd polar_to_complex(double r, double theta){
    Complexd result;
    result.re = r*cos(theta);
	result.im = r*sin(theta);
    return result;
}
Complexd add_complex(Complexd a, Complexd b){
    Complexd result;
    result.re = a.re + b.re;
    result.im = a.im + b.im;
    return result;
}
Complexd sub_complex(Complexd a, Complexd b){
	Complexd result;
	result.re = a.re - b.re;
	result.im = a.im - b.im;
	return result;
}
Complexd mult_complex(Complexd a, Complexd b){
	Complexd result;
	result.re = (a.re*b.re) - (a.im*b.im);
    result.im = (a.re*b.im) + (a.im*b.re);
	return result;
}

void fft_calc(int N,double *x,Complexd *X,Complexd *P,int step,Complexd *twids){
    Complexd *S = P + N/2;
    if (N == 1){
		X[0].re = x[0];
        X[0].im = 0.0;
		return;
    }
    
    fft_calc(N/2, x,      S,   X,2*step, twids);
    fft_calc(N/2, x+step, P,   X,2*step, twids);

    int k;
    for (k=0;k<N/2;k++){
		P[k] = mult_complex(P[k],twids[k*step]);
		X[k]     = add_complex(S[k],P[k]);
		X[k+N/2] = sub_complex(S[k],P[k]);
    }

}


int fft(double *x, int N, Complexd *X){

    Complexd *twiddle_factors = (Complexd*)malloc(sizeof(Complexd*)*(N/2));
    Complexd *Xt = (Complexd*)malloc(sizeof(Complexd)*N);
    int k;
    for (k=0;k<N/2;k++){
		twiddle_factors[k] = polar_to_complex(1.0, 2.0*PI*k/N);
    }
    fft_calc(N, x, X, Xt, 1, twiddle_factors);

    free(twiddle_factors);
    free(Xt);

    return 0;

}
