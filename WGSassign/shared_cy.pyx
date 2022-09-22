# cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True
import numpy as np
cimport numpy as np
from cython.parallel import prange
from cython import boundscheck, wraparound
from libc.math cimport sqrt

##### Shared Cython functions #####
# EM MAF update
cpdef emMAF_update(float[:,::1] L, float[::1] f, int t):
    cdef int m = L.shape[0]
    cdef int n = L.shape[1]//2
    cdef int i, s
    cdef float tmp, p0, p1, p2
    with nogil:
        for s in prange(m, num_threads=t):
            tmp = 0.0
            for i in range(n):
                p0 = L[s,2*i+0]*(1 - f[s])*(1 - f[s])
                p1 = L[s,2*i+1]*2*f[s]*(1 - f[s])
                p2 = (1.0 - L[s,2*i+0] - L[s,2*i+1])*f[s]*f[s]
                tmp = tmp + (p1 + 2*p2)/(2*(p0 + p1 + p2))
            f[s] = tmp/(<float>n)

# Root mean squared error (1D)
cpdef rmse1d(float[::1] v1, float[::1] v2):
    cdef int n = v1.shape[0]
    cdef int i
    cdef float res = 0.0
    for i in range(n):
        res = res + (v1[i] - v2[i])*(v1[i] - v2[i])
    res = res/(<float>n)
    return sqrt(res)
