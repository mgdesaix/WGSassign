# cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True
import numpy as np
cimport numpy as np
from cython.parallel import prange
from cython import boundscheck, wraparound

# Calculate zscore
# L is the reader_cy read beagle file, matrix M x (2 * N)
cpdef AD_summary(float[:,::1] AD, float[:,::1] L, int t, int i):
    cdef int m = L.shape[0]
    cdef int s
    cdef dict AD_GL_dict = {}
    cdef float like0, like1, like2
    with nogil:
            for s in prange(m, num_threads=t):
                key = tuple(AD[s,2*i], AD[s,2*i+1])
                if key not in AD_GL_dict.keys():
                    AD_GL_dict[key] = [[L[s,2*i], L[s,2*i+1], (1 - L[s,2*i] - L[s,2*i+1])]]
                else:
                    AD_GL_dict[key].append([L[s,2*i], L[s,2*i+1], (1 - L[s,2*i] - L[s,2*i+1])])
    return AD_GL_dict