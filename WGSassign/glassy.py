"""
Get likelihood of assignment
"""

__author__ = "Matt DeSaix"

# libraries
import numpy as np

# Import scripts
from WGSassign import glassy_cy

# log-likelihoods for assignment
# L is the GLs from beagle file (M x (N*2))
# A is the allele frequencies of the reference pops (M x K)
def assignLL(L, A, t):
    ## Initiate variables and containers
    # number of loci
    m = L.shape[0]
    # number of individuals
    n = L.shape[1] // 2
    # number of populations
    k = A.shape[1]
    # loglike matrix, of n x k (rows = individuals, columns = reference pops)
    logl_mat = np.zeros((n,k), dtype=np.float32)
    
    print(str(n) + " individuals to assign to " + str(k) + " populations")
    
    for i in range(n):
        for j in range(k):
            # set log like vector for individual i to pop j
            logl_vec = np.zeros(m, dtype=np.float32)
            # Check what's going on with double vs float
            print("A is " + str(np.array(list(A)).dtype))
            print("L is " + str(np.array(list(L)).dtype))
            print("logl_vec is " + str(np.array(list(logl_vec)).dtype))
            # fill vector
            glassy_cy.loglike(L, A, logl_vec, t, i, j)
            # loglike sum
            loglike = np.sum(logl_vec, dtype=float)
            print("Individual " + str(i) + " done for pop " + str(j))
            print("Log-likelihood: " + str(loglike))
            # fill output matrix
            logl_mat[i,j] = loglike
    del logl_vec
    return logl_mat

# EM algorithm for mixing proportions
# def emProp():
