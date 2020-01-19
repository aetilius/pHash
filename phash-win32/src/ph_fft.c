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


void fft_calc(const int N, const double *x, Complexd *X, Complexd *P, const int step, const Complexd *twids){
    int k;
    Complexd *S = P + N/2;
    if (N == 1){
		X[0].re = x[0];
        X[0].im = 0;
		return;
    }
    
    fft_calc(N/2, x,      S,   X,2*step, twids);
    fft_calc(N/2, x+step, P,   X,2*step, twids);
 
    for (k=0;k<N/2;k++){
		P[k] = mult_complex(P[k],twids[k*step]);
		X[k]     = add_complex(S[k],P[k]);
		X[k+N/2] = sub_complex(S[k],P[k]);
    }

}


int fft(const double *x, const int N, Complexd *X){

    Complexd *twiddle_factors = (Complexd*)malloc(sizeof(Complexd)*(N/2));
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
