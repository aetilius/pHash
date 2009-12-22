#include "phcomplex.h"

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

