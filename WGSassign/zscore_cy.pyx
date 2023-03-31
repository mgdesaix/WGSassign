# cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True
import numpy as np
cimport numpy as np
from cython.parallel import prange
from cython import boundscheck, wraparound
from libc.math cimport log

# Calculate zscore
# L is the reader_cy read beagle file, matrix M x (2 * N)
cpdef expected_W_l(float[:,::1] L, int[::1] L_keep, float[:,::1] A, long[:,::1] AD, int[:,::1] AD_array, float[:,::1] AD_factorial, float[:,::1] AD_like, int t, int i, int k, float[::1] W_l_obs_list, float[::1] W_l):
    cdef int m = L.shape[0]
    cdef int s, s_index, Dl, Aa, Ar
    cdef float A_sk, P_gl0, P_gl1, P_gl2, f_gl0, f_gl1, f_gl2, f_gl_log
    with nogil:
            for s_index in prange(m, num_threads=t):
                s = L_keep[s_index]
                A_sk = A[s,k]
                P_gl0 = (1-A_sk)*(1-A_sk)
                P_gl1 = 2*(1-A_sk)*A_sk
                P_gl2 = A_sk*A_sk
                f_gl0 = L[s,2*i+0] * P_gl0
                f_gl1 = L[s,2*i+1] * P_gl1
                f_gl2 = (1-L[s,2*i+0]-L[s,2*i+1]) * P_gl2
                f_gl_log = log(f_gl0 + f_gl1 + f_gl2)
                W_l_obs_list[s_index] = W_l_obs_list[s_index] + f_gl_log
                Dl = AD[s,2*i] + AD[s,2*i+1]
                for Aa in range(Dl+1):
                    Ar = Dl - Aa
                    ad_index = np.argwhere((AD_array[:,0] == Ar) & (AD_array[:,1] == Aa))[0][0]
                    W_l[s_index] = W_l[s_index] + f_gl_log * P_gl0 * AD_factorial[ad_index,0] * AD_like[ad_index,0]
                    W_l[s_index] = W_l[s_index] + f_gl_log * P_gl1 * AD_factorial[ad_index,0] * AD_like[ad_index,1]
                    W_l[s_index] = W_l[s_index] + f_gl_log * P_gl2 * AD_factorial[ad_index,0] * AD_like[ad_index,2]
                    
                    
                    
                    
                    