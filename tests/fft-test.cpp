#include "ph_fft.h"
#include <cassert>

using namespace std;

int main(int argc, char **argv){

    int N = 4;
    double signal[4] = { 1, 4, 3, 2};

	complex<double> X[N];
    int rc = fft(signal, N, X);
	assert(rc == 0);

    assert(abs(real(X[0]) - 10.0) < 0.0000001);
    assert(abs(imag(X[0]) - 0.0)  < 0.0000001);

    assert(abs(real(X[1]) + 2.0) < 0.0000001);
    assert(abs(imag(X[1]) - 2.0) < 0.0000001);

    assert(abs(real(X[2]) + 2.0) < 0.0000001);
    assert(abs(imag(X[2]) - 0.0) < 0.0000001);

    assert(abs(real(X[3]) + 2.0) < 0.0000001);
    assert(abs(imag(X[3]) + 2.0) < 0.0000001);

    double signal2[4] = { 3, 5, 9, 2};
   
    rc = fft(signal2, N, X);
	assert(rc == 0);
    
    assert(abs(real(X[0]) - 19.0) < 0.0000001);
    assert(abs(imag(X[0]) - 0.0)  < 0.0000001);

    assert(abs(real(X[1]) + 6.0) < 0.0000001);
    assert(abs(imag(X[1]) - 3.0) < 0.0000001);

    assert(abs(real(X[2]) - 5.0) < 0.0000001);
    assert(abs(imag(X[2]) - 0.0) < 0.0000001);

    assert(abs(real(X[3]) + 6.0) < 0.0000001);
    assert(abs(imag(X[3]) + 3.0) < 0.0000001);

    return 0;
}
