# cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True
import os
import numpy as np
cimport numpy as np
from cython import boundscheck, wraparound
from cython.parallel import prange
from libcpp.vector cimport vector
#from libcpp.iostrem cimport cout scientific
from libc.string cimport strtok, strdup
from libc.stdlib cimport atof
from libc.stdio cimport FILE, fclose, fopen, fprintf

DTYPE = np.float32
ctypedef np.float32_t DTYPE_t

cpdef tuple readBeagle(str beagle):
    cdef int c = 0
    cdef int i, m, n, s
    cdef bytes line_bytes
    cdef str line_str
    cdef char* line
    cdef char* token
    cdef char* delims = "\t \n"
    cdef vector[vector[float]] L
    cdef vector[float] L_ind

    # New: store site names and sample names
    cdef list site_names = []
    cdef list sample_names = []

    with os.popen("gunzip -c " + beagle) as f:
        # --- Parse header line ---
        line_bytes = str.encode(f.readline())
        line = line_bytes
        token = strtok(line, delims)  # 'marker'
        token = strtok(NULL, delims)  # 'allele1'
        token = strtok(NULL, delims)  # 'allele2'
        c = 0
        while True:  # this loop is taking care of the header
            token = strtok(NULL, delims)
            if token == NULL:
                break
            c += 1
            # every 3rd token corresponds to a sample name
            if c % 3 == 1:
                sample_names.append(<str>token.decode())

        n = c  # total number of GL columns
        n_inds = n // 3

        # --- Parse data lines ---
        for line_str in f:
            line_bytes = str.encode(line_str)
            line = line_bytes

            token = strtok(line, delims)
            site_names.append(<str>token.decode())  # CHR_POS

            token = strtok(NULL, delims)  # skip allele1
            token = strtok(NULL, delims)  # skip allele2

            for i in range(n):
                if (i + 1) % 3 == 0:
                    token = strtok(NULL, delims)  # skip GL2
                else:
                    L_ind.push_back(atof(strtok(NULL, delims)))
            L.push_back(L_ind)
            L_ind.clear()

    m = L.size()
    cdef np.ndarray[DTYPE_t, ndim=2] L_np = np.empty((m, 2 * n_inds), dtype=DTYPE)
    cdef float *L_ptr
    for s in range(m):
        L_ptr = &L[s][0]
        L_np[s] = np.asarray(<float[:(2 * n_inds)]> L_ptr)

    return (L_np, sample_names, site_names)

# Convert PLINK bed format to Beagle format
cpdef convertBed(float[:,::1] L, unsigned char[:,::1] G, int G_len, float e, \
                    int m, int n, int t):
    cdef signed char[4] recode = [0, 9, 1, 2]
    cdef unsigned char mask = 3
    cdef unsigned char byte, code
    cdef int i, s, b, bytepart
    with nogil:
        for s in prange(m, num_threads=t):
            i = 0
            for b in range(G_len):
                byte = G[s,b]
                for bytepart in range(4):
                    code = recode[byte & mask]
                    if code == 0:
                        L[s,2*i+0] = e*e
                        L[s,2*i+1] = 2*e*(1 - e)
                    elif code == 1:
                        L[s,2*i+0] = (1 - e)*e
                        L[s,2*i+1] = (1 - e)*(1 - e) + e*e
                    elif code == 2:
                        L[s,2*i+0] = (1 - e)*(1 - e)
                        L[s,2*i+1] = 2*e*(1 - e)
                    else:
                        L[s,2*i+0] = 0.333333
                        L[s,2*i+1] = 0.333333
                    byte = byte >> 2
                    i = i + 1
                    if i == n:
                        break

# Array filtering
cpdef filterArrays(float[:,::1] L, float[::1] f, unsigned char[::1] mask):
    cdef int m = L.shape[0]
    cdef int n = L.shape[1]
    cdef int c = 0
    cdef int i, s
    for s in range(m):
        if mask[s] == 1:
            for i in range(n):
                L[c,i] = L[s,i] # Genotype likelihoods
            f[c] = f[s] # Allele frequency
            c += 1
