# cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True
import numpy as np
cimport numpy as np
from cython.parallel import prange
from cython import boundscheck, wraparound
from libc.math cimport log

# Calculate fisher information
# L is the reader_cy read beagle file, matrix M x (2 * N)
# A is allele frequency matrix of pop k
# f_pop is the output vector of fisher information across all loci for a population
cpdef loglike(float[:,::1] L, float[:,::1] A, int i, int n, float[::1] f_pop):
    cdef int m = L.shape[0]
    cdef int s, t
    cdef float u, n1, n2, term
    cdef float term_sum, th, g0, g1, g2
    with nogil:
            for s in prange(m, num_threads=t):
                term_sum = 0
                th = A[s,i]
                for t in range(n):
                    g0 = L[s,2*t+0]
                    g1 = L[s,2*t+1]
                    g2 = 1.0 - g0 - g1
                    u = (g0*th*th) + (g1*2.0*th*(1.0 - th)) + (g2*(1.0 - th)*(1.0 - th))
                    n1 = 2.0*(g0 + g2 - (2.0*g1))
                    n2 = ((1.0 - th)*n1) + (2.0*(g1 - g0))
                    term = -1.0 * ((n1/u) - ((n2/u) * (n2/u)))
                    term_sum = term_sum + term 
                f_pop[s] = f_pop[s] + term_sum
