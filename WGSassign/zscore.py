"""
Calculate the z-score
"""

__author__ = "Matt DeSaix"

# libraries
import numpy as np

# Import scripts
# from WGSassign import glassy_cy
# from WGSassign import emMAF

def AD_summary(ad_ind0, L_ind0):
  AD_GL_dict = {}
  for i in np.arange(ad_ind0.shape[0]):
    # The key is the allele depth combination
    key = tuple(ad_ind0[i,])
    if key not in AD_GL_dict.keys():
      AD_GL_dict[key] = [[L_ind0[i,0], L_ind0[i,1], (1 - L_ind0[i,0] - L_ind0[i,1])]]
    else:
      AD_GL_dict[key].append([L_ind0[i,0], L_ind0[i,1], (1 - L_ind0[i,0] - L_ind0[i,1])])
  AD_summary_dict = {}
  for key, value in AD_GL_dict.items():
    AD_summary_dict[key] = [len(value), np.mean(AD_GL_dict[key], axis = 0)]
  return AD_GL_dict, AD_summary_dict

def get_L_keep(ad_ind0, L_ind0, AD_summary_dict, n_threshold):
  L_keep = np.empty(ad_ind0.shape[0], dtype = np.int32)
  # n_threshold = 1000
  loci_kept = 0
  loci_removed_N = 0
  loci_removed_GL = 0
  loci_zero = 0
  for i in np.arange(ad_ind0.shape[0]):
    key = tuple(ad_ind0[i,])
    if key == (0,0):
      L_keep[i] = 0
      loci_zero += 1
    else:
      max_id = np.argwhere(AD_summary_dict[key][1] == np.max(AD_summary_dict[key][1]))[0][0]
      L_ind_full = [L_ind0[i,0], L_ind0[i,1], 1 - L_ind0[i,0] - L_ind0[i,1]]
      if AD_summary_dict[key][0] < n_threshold:
        L_keep[i] = 0
        loci_removed_N += 1
      elif np.abs(AD_summary_dict[key][1][max_id] - L_ind_full[max_id]) > 0.01:
        L_keep[i] = 0
        loci_removed_GL += 1
      else:
        L_keep[i] = 1
        loci_kept += 1
  return L_keep, loci_kept

def get_expected_W_l(L_ind0, L_keep, ad_ind0, mafs_pop0, AD_summary_dict, e=0.01):
  W_l_obs = 0
  W_l = np.zeros(ad_ind0.shape[0], dtype = np.float32)
  for i in np.arange(ad_ind0.shape[0]):
    if L_keep[i] == 1:
      m_i = mafs_pop0[i]
      P_gl = [(1-m_i)*(1-m_i), 2*m_i*(1-m_i), m_i*m_i]
      f_gl = [L_ind0[i,0] * P_gl[0],L_ind0[i,1] * P_gl[1],(1-L_ind0[i,0]-L_ind0[i,1]) * P_gl[2]]
      f_gl_log = np.log(f_gl[0] + f_gl[1] + f_gl[2])
      W_l_obs += f_gl_log
      key = tuple(ad_ind0[i,])
      # Getting the total depth
      Dl = np.sum(list(key))
      for Aa in np.arange(Dl+1):
        Ar = Dl - Aa
        A = (Ar, Aa)
        iter_key = tuple([Ar, Aa])
        # Probability of read depths, combinatorial
        ad_factorial = np.math.factorial(Dl) / (np.math.factorial(Aa)*np.math.factorial(Ar))
        P_r_a = [ad_factorial*((1-e)**Ar)*(e**Aa), ad_factorial*((1/2)**Dl), ad_factorial*((1-e)**Aa)*(e**Ar)]
        for j in range(3):
          W_l[i] += f_gl_log * P_gl[j] * P_r_a[j] * 1 * AD_summary_dict[iter_key][1][j]
  return W_l_obs, W_l

def get_var_W_l(L_ind0, L_keep, ad_ind0, mafs_pop0, AD_summary_dict, W_l, e=0.01):
  var_W_l = np.zeros(ad_ind0.shape[0], dtype = np.float32)
  # e = 0.01
  for i in np.arange(ad_ind0.shape[0]):
    if L_keep[i] == 1:
      m_i = mafs_pop0[i]
      P_gl = [(1-m_i)*(1-m_i), 2*m_i*(1-m_i), m_i*m_i]
      f_gl = [L_ind0[i,0] * P_gl[0], L_ind0[i,1] * P_gl[1], (1-L_ind0[i,0]-L_ind0[i,1]) * P_gl[2]]
      f_gl_log = np.log(f_gl[0] + f_gl[1] + f_gl[2])
      key = tuple(ad_ind0[i,])
      # Getting the total depth
      Dl = np.sum(list(key))
      for Aa in np.arange(Dl+1):
        Ar = Dl - Aa
        A = (Ar, Aa)
        iter_key = tuple([Ar, Aa])
        # Probability of read depths, combinatorial
        ad_factorial = np.math.factorial(Dl) / (np.math.factorial(Aa)*np.math.factorial(Ar))
        P_r_a = [ad_factorial*((1-e)**Ar)*(e**Aa), ad_factorial*((1/2)**Dl), ad_factorial*((1-e)**Aa)*(e**Ar)]
        for j in range(3):
          var_W_l[i] += (W_l[i]-f_gl_log)**2 * P_gl[j] * P_r_a[j] * 1 * AD_summary_dict[iter_key][1][j]
  return var_W_l
