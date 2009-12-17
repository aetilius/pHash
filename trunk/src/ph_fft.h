#ifndef _FFT_H
#define _FFT_H

#define PI 3.1415926535897932

#include <math.h>
#include <complex.h>

static int nbadds;
static int nbmults;

int fft(double *x, int N, complex double *X);

#endif
