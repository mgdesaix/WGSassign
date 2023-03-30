# cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True
import numpy as np
cimport numpy as np
from cython.parallel import prange
from cython import boundscheck, wraparound
from libc.math cimport log

# Calculate zscore
# L is the reader_cy read beagle file, matrix M x (2 * N)
cpdef expected_W_l(float[:,::1] L, float[::1] L_keep, float[:,::1] A, float[:,::1] AD, dict AD_summary_dict, int t, int i, int k, float W_l_obs, float[::1] W_l, float e):
    cdef int m = L.shape[0]
    cdef int s, s_index, Dl, Aa, Ar
    cdef float A_sk, P_gl0, P_gl1, P_gl2, f_gl0, f_gl1, f_gl2, f_gl_log, ad_factorial, P_r_a0, P_r_a1, P_r_a2
    cdef str key, iter_key
    with nogil:
            for s_index in prange(m, num_threads=t):
                s = L_keep[s_index]
                A_sk = A[s,k]
                P_gl0 = (1-A_sk)*(1-A_sk)
                P_gl1 = 2*(1-A_sk)*A_sk
                P_gl2 = A_sk*A_sk
                f_gl0 = L[s,2*i] * P_gl0
                f_gl1 = L[s,2*i+1] * P_gl1
                f_gl2 = (1-L[s,2*i]-L[s,2*i+1]) * P_gl2
                f_gl_log = log(f_gl0 + f_gl1 + f_gl2)
                W_l_obs = W_l_obs + f_gl_log
                with gil:
                    key = str([AD[s,2*i], AD[s,2*i+1]])
                    Dl = AD[s,2*i] + AD[s,2*i+1]
                    for Aa in range(Dl+1):
                        Ar = Dl - Aa
                        iter_key = str([Ar, Aa])
                        ad_factorial = np.math.factorial(Dl) / (np.math.factorial(Aa)*np.math.factorial(Ar))
                        P_r_a0 = ad_factorial*((1.0-e)**Ar)*(e**Aa)
                        P_r_a1 = ad_factorial*((0.5)**Dl)
                        P_r_a2 = ad_factorial*((1.0-e)**Aa)*(e**Ar)
                        W_l[s_index] = W_l[s_index] + f_gl_log * P_gl0 * P_r_a0 * AD_summary_dict[iter_key][1][0]
                        W_l[s_index] = W_l[s_index] + f_gl_log * P_gl1 * P_r_a1 * AD_summary_dict[iter_key][1][1]
                        W_l[s_index] = W_l[s_index] + f_gl_log * P_gl2 * P_r_a2 * AD_summary_dict[iter_key][1][2]
                    
                    
                    
                    
                    