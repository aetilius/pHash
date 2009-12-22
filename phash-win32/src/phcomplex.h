#ifndef _PH_COMPLEX_H
#define _PH_COMPLEX_H

typedef struct phcomplex {
     double re;
     double im;
} Complexd;

Complexd polar_to_complex(double r, double theta);

Complexd add_complex(Complexd a, Complexd b);

Complexd sub_complex(Complexd a, Complexd b);

Complexd mult_complex(Complexd a, Complexd b);

#endif
