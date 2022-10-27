"""
Functions for mixture proportions
"""

__author__ = "Matt DeSaix"

# libraries
import numpy as np

def em_mix(L_mat, L_mat_index, iter):
  # number of reference populations
  n_source = L_mat.shape[1]
  # populations delineated for which to calculate mixing proportions
  harvest_pops = np.unique(L_mat_index[:,1])
  n_harvest = len(harvest_pops)
  
  em_out = np.empty((n_harvest, n_source), dtype = np.float32)
  for i in range(n_harvest):
    # subset likelihoods by harvest populations
    pop_index = np.argwhere(L_mat_index[:,1] == harvest_pops[i]).reshape(-1)
    L_mat_pop = np.ascontiguousarray(L_mat[pop_index,:])
    n_ind = L_mat_pop.shape[0]
    # create empty matrix for iterations
    pi_em_iters = np.empty((iter, n_source))
    # initialize equal mixture proportions
    pi_mat = np.diag(np.full(n_source,1))/n_source
    for j in range(iter):
      # update assignment probabilities
      L_pi_mat = np.matmul(np.exp(L_mat_pop), pi_mat)
      L_pi_mat_sum = L_pi_mat/L_pi_mat.sum(axis = 1, keepdims = True)
      # update mixture proportions
      pi_vec = L_pi_mat_sum.sum(axis = 0, keepdims = True) / n_ind
      pi_em_iters[j,:] = pi_vec
      pi_mat = np.diag(pi_vec.reshape(-1))
    # save off the last EM iteration
    em_out[i,:] = pi_em_iters[-1,:]
    del pi_em_iters
  em_mix_out = np.hstack((harvest_pops.reshape((n_harvest, 1)), em_out))
  return em_mix_out

def mcmc_mix(L_mat, L_mat_index, iter):
  # number of reference populations
  n_source = L_mat.shape[1]
  # populations delineated for which to calculate mixing proportions
  harvest_pops = np.unique(L_mat_index[:,1])
  n_harvest = len(harvest_pops)
  
  mcmc_out = np.empty((n_harvest, n_source), dtype = np.float32)
  for i in range(n_harvest):
    # subset likelihoods by harvest populations
    pop_index = np.argwhere(L_mat_index[:,1] == harvest_pops[i]).reshape(-1)
    L_mat_pop = np.ascontiguousarray(L_mat[pop_index,:])
    n_ind = L_mat_pop.shape[0]
    # create empty matrix for iterations
    pi_mcmc_iters = np.empty((iter, n_source))
    # initialize equal mixture proportions
    pi_mat = np.diag(np.full(n_source,1))/n_source
    # for probability distributions
    rng = np.random.default_rng()
    
    for j in range(iter):
      # update assignment probabilities
      L_pi_mat = np.matmul(np.exp(L_mat_pop), pi_mat)
      L_pi_mat_sum = L_pi_mat/L_pi_mat.sum(axis = 1, keepdims = True)
      # update mixture proportions
      ## multinomial
      vals = rng.multinomial(1, L_pi_mat_sum, n_ind)
      counts = vals.sum(axis = 0) + 0.001 # avoid error if 0
      ## dirichlet
      pi_vec = rng.dirichlet(counts, 1)
      pi_mcmc_iters[j,:] = pi_vec
      pi_mat = np.diag(pi_vec.reshape(-1))
    # save off the last EM iteration
    mcmc_out[i,:] = pi_mcmc_iters[-1,:]
    del pi_em_iters
  mcmc_mix_out = np.hstack((harvest_pops.reshape((n_harvest, 1)), mcmc_out))
  return mcmc_mix_out

  
