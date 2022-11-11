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
cpdef fisher_obs(float[:,::1] L, float[:,::1] A, int t, int i, int n, float[::1] f_pop):
    cdef int m = A.shape[0]
    cdef int s, r
    cdef float u, n1, n2, term
    cdef float term_sum, th, g0, g1, g2
    with nogil:
            for s in prange(m, num_threads=t):
                term_sum = 0
                th = A[s,i]
                for r in range(n):
                    g0 = L[s,2*r+0]
                    g1 = L[s,2*r+1]
                    g2 = 1.0 - g0 - g1
                    u = (g2*(1.0-th)*(1.0-th)) + (g1*2.0*th*(1.0 - th)) + (g0*th*th)
                    n1 = 2.0*(g2 + g0 - (2.0*g1))
                    n2 = (th*n1) + (2.0*(g1 - g2))
                    term = -1.0 * ((n1/u) - ((n0/u) * (n0/u)))
                    term_sum = term_sum + term 
                f_pop[s] = f_pop[s] + term_sum

cpdef ne_obs(float[::1] f_pop, float[:,::1] A, int t, int i, int n, float[::1] ne_pop):
    cdef int m = A.shape[0]
    cdef int s
    cdef float th, i_obs, n_tilde
    with nogil:
            for s in prange(m, num_threads=t):
                th = A[s,i]
                i_obs = f_pop[s]
                n_tilde = 0.5 * i_obs * th * (1 - th)
                ne_pop[s] = ne_pop[s] + n_tilde
