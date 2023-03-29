"""
Calculate the z-score
"""

__author__ = "Matt DeSaix"

# libraries
import numpy as np
# from WGSassign import zscore_cy

def AD_summary(L, AD, t, i):
  AD_GL_dict = {}
  for s in np.arange(AD.shape[0]):
    # The key is the allele depth combination
    key = str([AD[s,2*i], AD[s,2*i+1]])
    if key not in AD_GL_dict.keys():
      AD_GL_dict[key] = [[L[s,2*i], L[s,2*i+1], (1 - L[s,2*i] - L[s,2*i])]]
    else:
      AD_GL_dict[key].append([L[s,2*i], L[s,2*i+1], (1 - L[s,2*i] - L[s,2*i])])
  # AD_GL_dict = zscore_cy.AD_summary(AD, L, t, i)
  AD_summary_dict = {}
  for key, value in AD_GL_dict.items():
    AD_summary_dict[key] = [len(value), np.mean(AD_GL_dict[key], axis = 0)]
  return AD_GL_dict, AD_summary_dict

def get_L_keep(L, AD, AD_summary_dict, n_threshold, t, i):
  L_keep = np.empty(AD.shape[0], dtype = np.int32)
  # n_threshold = 1000
  loci_kept = 0
  loci_removed_N = 0
  loci_removed_GL = 0
  loci_zero = 0
  for s in np.arange(AD.shape[0]):
    key = str([AD[s,2*i], AD[s,2*i+1]])
    Dl = AD[s,2*i] + AD[s,2*i+1]
    if Dl == 0:
      L_keep[s] = 0
      loci_zero += 1
    else:
      max_id = np.argwhere(AD_summary_dict[key][1] == np.max(AD_summary_dict[key][1]))[0][0]
      L_ind_full = [L[s,2*i], L[s,2*i+1], (1 - L[s,2*i] - L[s,2*i])]
      if AD_summary_dict[key][0] < n_threshold:
        L_keep[s] = 0
        loci_removed_N += 1
      elif np.abs(AD_summary_dict[key][1][max_id] - L_ind_full[max_id]) > 0.01:
        L_keep[s] = 0
        loci_removed_GL += 1
      else:
        L_keep[s] = 1
        loci_kept += 1
  return L_keep, loci_kept

def get_expected_W_l(L, L_keep, A, AD, AD_summary_dict, t, i, k):
  W_l_obs = 0
  W_l = np.zeros(AD.shape[0], dtype = np.float32)
  e = 0.01
  for s in np.arange(AD.shape[0]):
    if L_keep[s] == 1:
      A_sk = A[s,k]
      P_gl = [(1-A_sk)*(1-A_sk), 2*A_sk*(1-A_sk), A_sk*A_sk]
      f_gl = [L[s,2*i] * P_gl[0],L[s,2*i+1] * P_gl[1],(1-L[s,2*i]-L[s,2*i+1]) * P_gl[2]]
      f_gl_log = np.log(f_gl[0] + f_gl[1] + f_gl[2])
      W_l_obs += f_gl_log
      key = str([AD[s,2*i], AD[s,2*i+1]])
      # Getting the total depth
      Dl = AD[s,2*i] + AD[s,2*i+1]
      for Aa in np.arange(Dl+1):
        Ar = Dl - Aa
        # A = (Ar, Aa)
        iter_key = str([Ar, Aa])
        # iter_key = tuple([Ar, Aa])
        # Probability of read depths, combinatorial
        ad_factorial = np.math.factorial(Dl) / (np.math.factorial(Aa)*np.math.factorial(Ar))
        P_r_a = [ad_factorial*((1-e)**Ar)*(e**Aa), ad_factorial*((1/2)**Dl), ad_factorial*((1-e)**Aa)*(e**Ar)]
        for j in range(3):
          W_l[s] += f_gl_log * P_gl[j] * P_r_a[j] * 1 * AD_summary_dict[iter_key][1][j]
  return W_l_obs, W_l

def get_var_W_l(L, L_keep, A, AD, AD_summary_dict, W_l, t, i, k):
  var_W_l = np.zeros(AD.shape[0], dtype = np.float32)
  e = 0.01
  for s in np.arange(AD.shape[0]):
    if L_keep[s] == 1:
      A_sk = A[s,k]
      P_gl = [(1-A_sk)*(1-A_sk), 2*A_sk*(1-A_sk), A_sk*A_sk]
      f_gl = [L[s,2*i] * P_gl[0],L[s,2*i+1] * P_gl[1],(1-L[s,2*i]-L[s,2*i+1]) * P_gl[2]]
      f_gl_log = np.log(f_gl[0] + f_gl[1] + f_gl[2])
      key = str([AD[s,2*i], AD[s,2*i+1]])
      # Getting the total depth
      Dl = AD[s,2*i] + AD[s,2*i+1]
      for Aa in np.arange(Dl+1):
        Ar = Dl - Aa
        # A = (Ar, Aa)
        iter_key = str([Ar, Aa])
        # iter_key = tuple([Ar, Aa])
        # Probability of read depths, combinatorial
        ad_factorial = np.math.factorial(Dl) / (np.math.factorial(Aa)*np.math.factorial(Ar))
        P_r_a = [ad_factorial*((1-e)**Ar)*(e**Aa), ad_factorial*((1/2)**Dl), ad_factorial*((1-e)**Aa)*(e**Ar)]
        for j in range(3):
          var_W_l[s] += (W_l[s]-f_gl_log)**2 * P_gl[j] * P_r_a[j] * 1 * AD_summary_dict[iter_key][1][j]
  return var_W_l
