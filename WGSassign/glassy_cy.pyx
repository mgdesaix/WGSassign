# cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True
import numpy as np
cimport numpy as np
from cython.parallel import prange
from cython import boundscheck, wraparound
from libc.math cimport log

# Calculate assignment likelihoods
# L is the reader_cy read beagle file, matrix M x (2 * N)
# A is allele frequency matrix of pop k
# loglike_vec is the output vector of log likelihoods of assignment to the pop
cpdef loglike(float[:,::1] L, float[:,::1] A, float[::1] loglike_vec, int t, int i, int k):
    cdef int m = L.shape[0]
    cdef int s
    cdef float like0, like1, like2
    with nogil:
            for s in prange(m, num_threads=t):
                like0 = L[s,2*i+0]*(1 - A[s,k])*(1 - A[s,k])
                like1 = L[s,2*i+1]*2*(1 - A[s,k])*A[s,k]
                like2 = (1 - L[s,2*i+0] - L[s,2*i+1]) * A[s,k] * A[s,k]
                loglike_vec[s] = loglike_vec[s] + log(like0 + like1 + like2)
