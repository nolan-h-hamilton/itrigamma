# cython: language_level=3
import numpy as np
cimport numpy as np
from libc.stddef cimport size_t

cdef extern from "itrigamma.h":
    double itrigamma(double y)
    void itrigamma_vec(const double *y, double *out, size_t n)

cpdef double itrigamma_scalar(double y):
    return itrigamma(y)

def itrigamma_np(y):
    cdef np.ndarray[np.float64_t, ndim=1, mode="c"] y1
    cdef np.ndarray[np.float64_t, ndim=1, mode="c"] out1
    cdef size_t n
    cdef tuple shape

    arr = np.asarray(y, dtype=np.float64)
    shape = arr.shape
    y1 = np.ascontiguousarray(arr.ravel(), dtype=np.float64)
    out1 = np.empty_like(y1)

    n = <size_t> y1.size
    itrigamma_vec(<const double*> y1.data, <double*> out1.data, n)

    return out1.reshape(shape)
