#include <stdlib.h>
#include <math.h>

#ifndef _PH_COMPLEX_H
#define _PH_COMPLEX_H

typedef struct phcomplex {
     double re;
     double im;
} Complexd;

Complexd polar_to_complex(const double r, const double theta);

Complexd add_complex(const Complexd a, const Complexd b);

Complexd sub_complex(const Complexd a, const Complexd b);

Complexd mult_complex(const Complexd a, const Complexd b);

#endif
