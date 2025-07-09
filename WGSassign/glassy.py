"""
Get likelihood of assignment
"""

__author__ = "Matt DeSaix"

# libraries
import numpy as np

# Import scripts
from WGSassign import glassy_cy
from WGSassign import emMAF

# log-likelihoods for assignment
# L is the GLs from beagle file (M x (N*2))
# af is the allele frequencies of the reference pops (M x K)
def assignLL(L, af, t):
    ## Initiate variables and containers
    # number of loci
    m = L.shape[0]
    # number of individuals
    n = L.shape[1] // 2
    # number of populations
    k = af.shape[1]
    # loglike matrix, of n x k (rows = individuals, columns = reference pops)
    logl_mat = np.zeros((n,k), dtype=np.float32)
    
    print(str(n) + " individuals to assign to " + str(k) + " populations")
    
    for i in range(n):
        for j in range(k):
            # set log like vector for individual i to pop j
            logl_vec = np.zeros(m, dtype=np.float32)
            # fill vector
            glassy_cy.loglike(L, af, logl_vec, t, i, j)
            # loglike sum
            loglike = np.sum(logl_vec, dtype=float)
            # print("Individual " + str(i) + " done for pop " + str(j))
            # print("Log-likelihood: " + str(loglike))
            # fill output matrix
            logl_mat[i,j] = loglike
    del logl_vec
    return logl_mat

# Leave-one-out likelihoods
def loo(L, af, IDs, t, maf_iter, maf_tole, downsampled_L=None):
    ## Initiate variables and containers
    # number of loci
    m = L.shape[0]
    # number of individuals
    n = L.shape[1] // 2
    # number of populations
    k = af.shape[1]
    # loglike matrix, of n x k (rows = individuals, columns = reference pops)
    logl_mat = np.zeros((n,k), dtype=np.float32)
    
    print(str(n) + " individuals to assign to " + str(k) + " populations")
    if downsampled_L is not None:
        print("Using downsampled GLs for likelihood evaluation in LOO assignment.")

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
        af_pop = emMAF.emMAF(L_pop, maf_iter, maf_tole, t)
        # set minimum value for allele frequencies as 1 + the number of individuals sampled
        n_pop = L_pop.shape[1] // 2
        min_val = 1 / (2 * (n_pop + 1))
        max_val = 1 - min_val
        
        af_pop[af_pop < min_val] = min_val
        af_pop[af_pop > max_val] = max_val
        # column index of the allele frequency file (based on unique order of pops)
        pop_col = np.argwhere(pops == i_pop)[0,0]
        # put new allele frequencies in column
        af[:,pop_col] = af_pop
        del L_pop, af_pop
        
        for j in range(k):
            # set log like vector for individual i to pop j
            logl_vec = np.zeros(m, dtype=np.float32)
            # choose which GLs to use (for downsampling test option)
            L_source = downsampled_L if downsampled_L is not None else L
            # fill vector
            glassy_cy.loglike(L_source, af, logl_vec, t, i, j)
            # loglike sum
            loglike = np.sum(logl_vec, dtype=float)
            # print("Individual " + str(i) + " done for pop " + str(j))
            # print("Log-likelihood: " + str(loglike))
            # fill output matrix
            logl_mat[i,j] = loglike
    del logl_vec
    return logl_mat

