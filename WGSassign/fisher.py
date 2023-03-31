"""
Functions for Fisher information
"""

__author__ = "Matt DeSaix"

# libraries
import numpy as np
# from WGSassign import fisher_cy

def fisher_obs(L, af, IDs, t):
  # Number of loci
  m = L.shape[0]
  # Unique reference pop names
  pops = np.unique(IDs[:,1])
  # number of reference pops
  npops = len(pops)
  # create observed Fisher matrix
  f_obs = np.empty((m, npops), dtype=np.float32)
  # create effective sample size matrix
  ne_obs = np.empty((m, npops), dtype = np.float32)
  for i in range(npops):
    # get indices of which rows in ID file correspond to the given reference pop
    pop_index = np.argwhere(IDs[:,1] == pops[i])
    # convert indices to relevant column indices of "L" file in pcangsd (Remember Beagle file is converted to 2 cols per individual)
    L1 = pop_index*2
    L2 = L1 + 1
    L_cat = np.concatenate((L1, L2))
    L_cat_index = np.sort(L_cat, axis = 0).reshape(-1)
    L_pop = np.ascontiguousarray(L[:,L_cat_index])
    # specify population allele frequency column from matrix
    # af_pop = af[:,i]
    f_pop = np.zeros(m, dtype=np.float32)
    # calculate fisher information for given reference pop i with n individuals
    n = L_pop.shape[1]//2
    fisher_cy.fisher_obs(L_pop, af, t, i, n, f_pop)
    f_obs[:,i] = f_pop
    # calculate effective sample size from fisher information
    ne_pop = np.zeros(m, dtype=np.float32)
    fisher_cy.ne_obs(f_pop, af, t, i, n, ne_pop)
    ne_obs[:,i] = ne_pop
    del L_pop, f_pop, ne_pop
  return f_obs, ne_obs

def fisher_obs_ind(L, af, IDs, t):
  # Number of loci
  m = L.shape[0]
  # Number of individuals
  n = L.shape[1] // 2
  pops = np.unique(IDs[:,1])
  ne_ind_full = np.zeros(n, dtype=np.float32)
  for i in range(n):
    f_ind = np.zeros(m, dtype=np.float32)
    pop_i = np.argwhere(pops == IDs[i,1])[0][0]
    fisher_cy.fisher_obs_ind(L, af, t, i, pop_i, f_ind)
    ne_ind = np.zeros(m, dtype=np.float32)
    fisher_cy.ne_obs_ind(f_ind, af, t, pop_i, ne_ind)
    ne_ind_full[i] = ne_ind_full[i] + np.mean(ne_ind)
  return ne_ind_full
  
