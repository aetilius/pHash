# distutils: language=c++
import os
from libc.stdlib cimport free
from libc.stdint cimport uint8_t, uint64_t

from phash cimport (ph_mh_imagehash, ph_hammingdistance2,
                    ph_dct_imagehash, ph_hamming_distance)


cdef class MHImageHash:
    cdef bytes c_hash

    def __getstate__(self):
        return self.c_hash

    def __setstate__(self, state):
        self.c_hash = state

    def hamming_distance(self, MHImageHash other):
        return int(round(100 * ph_hammingdistance2(self.c_hash, len(self.c_hash), other.c_hash, len(other.c_hash))))

    @staticmethod
    def from_path(path, alpha=2.0, level=1.0):
        cdef uint8_t* c_hash
        cdef int c_hashlen = 0
        c_hash = ph_mh_imagehash(os.fsencode(path), c_hashlen, alpha, level)
        if c_hash is NULL:
            raise RuntimeError("ph_mh_imagehash failed")
        else:
            obj = MHImageHash()
            try:
                obj.c_hash = c_hash[:c_hashlen]
            finally:
                free(c_hash)
            return obj


cdef class DCTImageHash:
    cdef uint64_t c_hash

    def __getstate__(self):
        return self.c_hash

    def __setstate__(self, state):
        self.c_hash = state

    def hamming_distance(self, DCTImageHash other):
        return ph_hamming_distance(self.c_hash, other.c_hash)

    @staticmethod
    def from_path(path):
        obj = DCTImageHash()
        if ph_dct_imagehash(os.fsencode(path), obj.c_hash) == 0:
            return obj
        else:
            raise RuntimeError("ph_dct_imagehash failed")

