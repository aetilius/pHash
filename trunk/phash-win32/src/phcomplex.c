#include "phcomplex.h"

Complexd polar_to_complex(const double r, const double theta){
    Complexd result;
    result.re = r*cos(theta);
	result.im = r*sin(theta);
    return result;
}
Complexd add_complex(const Complexd a, const Complexd b){
    Complexd result;
    result.re = a.re + b.re;
    result.im = a.im + b.im;
    return result;
}
Complexd sub_complex(const Complexd a, const Complexd b){
	Complexd result;
	result.re = a.re - b.re;
	result.im = a.im - b.im;
	return result;
}
Complexd mult_complex(const Complexd a, const Complexd b){
	Complexd result;
	result.re = (a.re*b.re) - (a.im*b.im);
    	result.im = (a.re*b.im) + (a.im*b.re);
	return result;
}

