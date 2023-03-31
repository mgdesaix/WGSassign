"""
Calculate the z-score
"""

__author__ = "Matt DeSaix"

# libraries
import numpy as np
# from WGSassign import zscore_cy

def AD_summary(L, AD, i, n_threshold):
  AD_GL_dict = {}
  for s in np.arange(AD.shape[0]):
    # The key is the allele depth combination
    key = tuple([AD[s,2*i], AD[s,2*i+1]])
    if key not in AD_GL_dict.keys():
      AD_GL_dict[key] = [[L[s,2*i], L[s,2*i+1], (1 - L[s,2*i] - L[s,2*i])]]
    else:
      AD_GL_dict[key].append([L[s,2*i], L[s,2*i+1], (1 - L[s,2*i] - L[s,2*i])])
  AD_summary_dict = {}
  for key, value in AD_GL_dict.items():
    AD_summary_dict[key] = [len(value), np.mean(AD_GL_dict[key], axis = 0)]
  AD_summary_array = np.empty((len(AD_summary_dict.keys()), 4), dtype = np.int32)
  for j in range(len(list(AD_summary_dict.keys()))):
    key = list(AD_summary_dict.keys())[j]
    a1 = list(list(AD_summary_dict.keys())[j])[0]
    a2 = list(list(AD_summary_dict.keys())[j])[1]
    n_loci = AD_summary_dict[key][0]
    AD_summary_array[j,:] = [a1, a2, a1+a2, n_loci]
  AD_filtered = AD_summary_array[(AD_summary_array[:,3] > n_threshold) & (AD_summary_array[:,2] != 0)]
  dict_sum = AD_filtered[:,0] + AD_filtered[:,1]
  dl, dl_counts = np.unique(dict_sum, return_counts=True)
  dl_keep = dl[dl < dl_counts]
  AD_array = AD_filtered[np.in1d(AD_filtered[:,2], dl_keep)]
    
  return AD_summary_dict, AD_array

def get_L_keep(L, AD, AD_summary_dict, AD_array, i):
  L_keep = np.empty(AD.shape[0], dtype = np.int32)
  for s in np.arange(AD.shape[0]):
    Ar = AD[s,2*i]
    Aa = AD[s,2*i+1]
    observed_A = len(np.argwhere((AD_array[:,0] == Ar) & (AD_array[:,1] == Aa)))
    if observed_A != 1:
      L_keep[s] = 0
    else:
      key = tuple([AD[s,2*i], AD[s,2*i+1]])
      max_id = np.argwhere(AD_summary_dict[key][1] == np.max(AD_summary_dict[key][1]))[0][0]
      L_ind_full = [L[s,2*i], L[s,2*i+1], (1 - L[s,2*i] - L[s,2*i])]
      if np.abs(AD_summary_dict[key][1][max_id] - L_ind_full[max_id]) > 0.01:
        L_keep[s] = 0
      else:
        L_keep[s] = 1
  L_keep_final = np.argwhere(L_keep == 1).reshape(-1).astype(np.int32)
  loci_kept = L_keep_final.shape[0]
  return L_keep_final, loci_kept

def get_factorials(AD_array, AD_summary_dict, e):
  AD_factorial = np.zeros((AD_array.shape[0], 3), dtype = float)
  AD_like = np.zeros((AD_array.shape[0], 3), dtype = float)
  for i in range(AD_factorial.shape[0]):
    Ar = AD_array[i,0]
    Aa = AD_array[i,1]
    Dl = Aa + Ar
    ad_factorial = np.math.factorial(Dl) / (np.math.factorial(Aa)*np.math.factorial(Ar))
    P_r_a0 = ad_factorial*((1-e)**Ar)*(e**Aa)
    P_r_a1 = ad_factorial*((1/2)**Dl)
    P_r_a2 = ad_factorial*((1-e)**Aa)*(e**Ar)
    AD_factorial[i,:] = [P_r_a0, P_r_a1, P_r_a2]
    AD_like[i:] = AD_summary_dict[tuple([Ar, Aa])][1]
  return AD_factorial, AD_like

def get_expected_W_l(L, L_keep, A, AD, AD_array, AD_factorial, AD_like, t, i, k):
  W_l_obs_list = np.zeros(L_keep.shape[0], dtype = np.float32)
  W_l = np.zeros(L_keep.shape[0], dtype = np.float32)
  # zscore_cy.expected_W_l(L, L_keep, A, AD, AD_summary_dict, t, i, k, e, W_l_obs_list, W_l)
  for s_index in range(L_keep.shape[0]):
    s = L_keep[s_index]
    A_sk = A[s,k]
    P_gl = [(1-A_sk)*(1-A_sk), 2*A_sk*(1-A_sk), A_sk*A_sk]
    f_gl = [L[s,2*i] * P_gl[0],L[s,2*i+1] * P_gl[1],(1-L[s,2*i]-L[s,2*i+1]) * P_gl[2]]
    f_gl_log = np.log(f_gl[0] + f_gl[1] + f_gl[2])
    W_l_obs_list[s_index] = f_gl_log
    # Getting the total depth
    Dl = AD[s,2*i] + AD[s,2*i+1]
    for Aa in np.arange(Dl+1):
      Ar = Dl - Aa
      ad_index = np.argwhere((AD_summary_array_filtered[:,0] == Ar) & (AD_summary_array_filtered[:,1] == Aa))[0][0]
      W_l[s_index] = W_l[s_index] + f_gl_log * P_gl[0] * AD_factorial[ad_index,0] * AD_like[ad_index,0]
      W_l[s_index] = W_l[s_index] + f_gl_log * P_gl[1] * AD_factorial[ad_index,1] * AD_like[ad_index,1]
      W_l[s_index] = W_l[s_index] + f_gl_log * P_gl[2] * AD_factorial[ad_index,2] * AD_like[ad_index,2]
  W_l_obs = np.sum(W_l_obs_list, dtype=float)
  return W_l_obs, W_l

def get_var_W_l(L, L_keep, A, AD, AD_array, AD_factorial, AD_like, W_l, t, i, k):
  var_W_l = np.zeros(L_keep.shape[0], dtype = np.float32)
  for s_index in range(L_keep.shape[0]):
    s = L_keep[s_index]
    A_sk = A[s,k]
    P_gl = [(1-A_sk)*(1-A_sk), 2*A_sk*(1-A_sk), A_sk*A_sk]
    f_gl = [L[s,2*i] * P_gl[0],L[s,2*i+1] * P_gl[1],(1-L[s,2*i]-L[s,2*i+1]) * P_gl[2]]
    f_gl_log = np.log(f_gl[0] + f_gl[1] + f_gl[2])
    # Getting the total depth
    Dl = AD[s,2*i] + AD[s,2*i+1]
    for Aa in np.arange(Dl+1):
      Ar = Dl - Aa
      ad_index = np.argwhere((AD_summary_array_filtered[:,0] == Ar) & (AD_summary_array_filtered[:,1] == Aa))[0][0]
      var_W_l[s_index] = var_W_l[s_index] + (W_l[s_index]-f_gl_log)**2 * P_gl[0] * AD_factorial[ad_index,0] * AD_like[ad_index,0]
      var_W_l[s_index] = var_W_l[s_index] + (W_l[s_index]-f_gl_log)**2 * P_gl[1] * AD_factorial[ad_index,1] * AD_like[ad_index,1]
      var_W_l[s_index] = var_W_l[s_index] + (W_l[s_index]-f_gl_log)**2 * P_gl[2] * AD_factorial[ad_index,2] * AD_like[ad_index,2]
  return var_W_l
