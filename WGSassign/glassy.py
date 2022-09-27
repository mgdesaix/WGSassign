"""
Get likelihood of assignment
"""

__author__ = "Matt DeSaix"

# libraries
import numpy as np

# Import scripts
from WGSassign import glassy_cy
from WGSassign import shared

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
            # fill vector
            glassy_cy.loglike(L, A, logl_vec, t, i, j)
            # loglike sum
            loglike = np.sum(logl_vec, dtype=float)
            # print("Individual " + str(i) + " done for pop " + str(j))
            # print("Log-likelihood: " + str(loglike))
            # fill output matrix
            logl_mat[i,j] = loglike
    del logl_vec
    return logl_mat

# Leave-one-out likelihoods
def loo(L, A, IDs, t, maf_iter, maf_tole):
    ## Initiate variables and containers
    # number of loci
    m = L.shape[0]
    # number of individuals
    n = L.shape[1] // 2
    # number of populations
    k = A.shape[1]
    # loglike matrix, of n x k (rows = individuals, columns = reference pops)
    logl_mat = np.zeros((n,k), dtype=np.float32)
    # set minimum value for allele frequencies as 1 + the number of individuals sampled
    # min_val = 1 / (2 * (n + 1))
    min_val = 1e-10
    
    print(str(n) + " individuals to assign to " + str(k) + " populations")
    # unique reference pops
    pops = np.unique(IDs[:,1])
    for i in range(n):
        # for each individual, recalculate allele frequency
        i_pop = IDs[i,1]
        # indices of all the individuals in the pop
        pop_index = np.argwhere(IDs[:,1] == i_pop)
        # remove current individual i
        pop_index_except = pop_index[pop_index != i]
        
        L1 = pop_index_except * 2
        L2 = L1 + 1
        L_cat = np.concatenate((L1, L2))
        L_cat_index = np.sort(L_cat, axis = 0).reshape(-1)
        L_pop = np.ascontiguousarray(L[:,L_cat_index])
        f_pop = shared.emMAF(L_pop, maf_iter, maf_tole, t)
        f_pop[f_pop < min_val] = min_val
        # column index of the allele frequency file (based on unique order of pops)
        pop_col = np.argwhere(pops == i_pop)[0,0]
        # put new allele frequencies in column
        A[:,pop_col] = f_pop
        del L_pop, f_pop
        
        for j in range(k):
            # set log like vector for individual i to pop j
            logl_vec = np.zeros(m, dtype=np.float32)
            # fill vector
            glassy_cy.loglike(L, A, logl_vec, t, i, j)
            # loglike sum
            loglike = np.sum(logl_vec, dtype=float)
            # print("Individual " + str(i) + " done for pop " + str(j))
            # print("Log-likelihood: " + str(loglike))
            # fill output matrix
            logl_mat[i,j] = loglike
    del logl_vec
    return logl_mat

# EM algorithm for mixing proportions
# def emProp():
