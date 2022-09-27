# distutils: language=c++
from libc.stdint cimport uint8_t, uint64_t

cdef extern from "../../src/pHash.h":
    uint8_t* ph_mh_imagehash(const char *, int&, float, float)
    double ph_hammingdistance2(uint8_t*, int, uint8_t*, int)
    int ph_dct_imagehash(const char *path, uint64_t&)
    int ph_hamming_distance(uint64_t hash1, uint64_t hash1)
