"""
WGSassign.
Main caller.

Matt DeSaix
"""

__author__ = "Matt DeSaix"

# Libraries
import argparse
import os
import subprocess
import sys
from datetime import datetime

# Find length of PLINK files
def extract_length(filename):
	process = subprocess.Popen(['wc', '-l', filename], stdout=subprocess.PIPE)
	result, err = process.communicate()
	return int(result.split()[0])

# Argparse
parser = argparse.ArgumentParser(prog="WGSassign")
parser.add_argument("-b", "--beagle", metavar="FILE",
	help="Filepath to genotype likelihoods in gzipped Beagle format from ANGSD")
# parser.add_argument("-p", "--plink", metavar="FILE-PREFIX",
# 	help="Prefix PLINK files (.bed, .bim, .fam)")
parser.add_argument("-t", "--threads", metavar="INT", type=int, default=1,
	help="Number of threads")
parser.add_argument("-o", "--out", metavar="OUTPUT", default="wgsassign",
	help="Prefix for output files")
parser.add_argument("--maf_iter", metavar="INT", type=int, default=200,
	help="Maximum iterations for minor allele frequencies estimation - EM (200)")
parser.add_argument("--maf_tole", metavar="FLOAT", type=float, default=1e-4,
	help="Tolerance for minor allele frequencies estimation update - EM (1e-4)")
# parser.add_argument("--iter", metavar="INT", type=int, default=100,
# 	help="Maximum iterations for estimation of individual allele frequencies (100)")
# parser.add_argument("--tole", metavar="FLOAT", type=float, default=1e-5,
# 	help="Tolerance for update in estimation of individual allele frequencies (1e-5)")

#############################################################3
# Reference population allele frequencies
parser.add_argument("--pop_af_IDs", metavar="FILE",
	help="Filepath to IDs for reference population beagle")
parser.add_argument("--get_reference_af", action="store_true", 
  help="Estimate allele frequencies for reference populations")

# Individual effective sample size
parser.add_argument("--ne_obs_ind", action="store_true", 
  help="Estimate individuals effective sample sizes")

# Estimate likelihoods of assignment
parser.add_argument("--pop_af_file", metavar="FILE",
	help="Filepath to reference population allele frequencies")
parser.add_argument("--get_pop_like", action="store_true", 
  help="Estimate log likelihood of individual assignment to each reference population")

# Leave one out assignment accuracy
parser.add_argument("--loo", action="store_true",
	help="Perform leave-one-out cross validation")

# z-score
parser.add_argument("--get_assignment_z_score", action="store_true", 
  help="Calculate z-score for individuals")
parser.add_argument("--get_reference_z_score", action="store_true", 
  help="Calculate z-score for individuals")
parser.add_argument("--ind_ad_file", metavar="FILE",
	help="Filepath to individual allele depths")
parser.add_argument("--allele_count_threshold", metavar="INT", type=int,
	help="Minimum number of loci needed to keep a specific allele count combination")
parser.add_argument("--single_read_threshold", action="store_true",
	help="Use only loci with a single read. Helpful for computational efficiency when individuals's sequencing depths vary.")
parser.add_argument("--ind_start", metavar="INT", type=int,
	help="Start analysis at this individual index (0-index: i.e. 0 starts with the 1st individual)")
parser.add_argument("--ind_end", metavar="INT", type=int,
	help="End analysis at this individual index (0-index: i.e. If you have 10 individuals, 9 is the 10th individual)")

# Mixture proportions
parser.add_argument("--pop_like", metavar="FILE",
	help="Filepath to population assignment log likelihood file")
parser.add_argument("--pop_like_IDs", metavar="FILE",
	help="Filepath to IDs for population assignment log likelihood file")
parser.add_argument("--get_em_mix", action="store_true", 
  help="Estimate mixture proportions with EM algorithm")
parser.add_argument("--get_mcmc_mix", action="store_true", 
  help="Estimate mixture proportions with MCMC algorithm")
parser.add_argument("--mixture_iter", metavar="INT", type=int, default=200,
	help="Maximum iterations mixture estimation - EM (200)")

###################################################################

##### WGSassign #####
def main():
	args = parser.parse_args()
	if len(sys.argv) < 2:
		parser.print_help()
		sys.exit()
	print("WGSassign")
	print("Matt DeSaix.")
	print("Using " + str(args.threads) + " thread(s).\n")

	# Check input
	# assert (args.beagle is not None) or (args.plink is not None), \
	# 		"Please provide input data (args.beagle or args.plink)!"

	# Create log-file of arguments
	full = vars(parser.parse_args())
	deaf = vars(parser.parse_args([]))
	with open(args.out + ".args", "w") as f:
		f.write("WGSassign\n")
		f.write("Time: " + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + "\n")
		f.write("Directory: " + str(os.getcwd()) + "\n")
		f.write("Options:\n")
		for key in full:
			if full[key] != deaf[key]:
				if type(full[key]) is bool:
					f.write("\t-" + str(key) + "\n")
				else:
					f.write("\t-" + str(key) + " " + str(full[key]) + "\n")
	del full, deaf

	# Control threads
	os.environ["OMP_NUM_THREADS"] = str(args.threads)
	os.environ["OPENBLAS_NUM_THREADS"] = str(args.threads)
	os.environ["MKL_NUM_THREADS"] = str(args.threads)

	# Numerical libraries
	import numpy as np
	from math import ceil

	# Import scripts
	from WGSassign import emMAF
	from WGSassign import glassy
	from WGSassign import mixture
	from WGSassign import fisher
	from WGSassign import reader_cy
	from WGSassign import zscore

	# Parse data
	if args.beagle is not None:
		print("Parsing Beagle file.")
		assert os.path.isfile(args.beagle), "Beagle file doesn't exist!"
		L = reader_cy.readBeagle(args.beagle)
		m = L.shape[0]
		n = L.shape[1]//2
		print("Loaded " + str(m) + " sites and " + str(n) + " individuals.")
	
####################################
  # Reference population allele frequencies
	if args.get_reference_af:
	  print("Parsing reference population ID file.")
	  assert os.path.isfile(args.pop_af_IDs), "Reference population ID file does not exist!!"
	  # File is tab-delimited with 2 columns and each row is an individual
	  # 1st column = Individual sample names corresponding to Beagle file
	  # 2nd column = Reference pop name
	  IDs = np.loadtxt(args.pop_af_IDs, delimiter = "\t", dtype = "str")
	  # Unique reference pop names
	  pops = np.unique(IDs[:,1])
	  # number of reference pops
	  npops = len(pops)
	  m = L.shape[0] # number of sites
	  # allele frequency matrix
	  af = np.empty((m, npops), dtype=np.float32)
	  # number of individuals from beagle
	  n = L.shape[1] // 2
	  # Check number of individuals from beagle is same as reference file
	  assert (n == IDs.shape[0]), "Number of individuals in beagle and reference ID file do not match!"
	  
	  # For each reference population, estimate the allele frequencies from the beagle file
	  for i in range(npops):
	    # get indices of which rows in ID file correspond to the given reference pop
	    pop_index = np.argwhere(IDs[:,1] == pops[i])
	    # convert indices to relevant column indices of "L" file in pcangsd (Remember Beagle file is converted to 2 cols per individual)
	    L1 = pop_index*2
	    L2 = L1 + 1
	    L_cat = np.concatenate((L1, L2))
	    L_cat_index = np.sort(L_cat, axis = 0).reshape(-1)
	    L_pop = np.ascontiguousarray(L[:,L_cat_index])
	    af_pop = emMAF.emMAF(L_pop, args.maf_iter, args.maf_tole, args.threads)
	    # set minimum value for allele frequencies as 1 + the number of individuals sampled
	    n_pop = L_pop.shape[1] // 2
	    min_val = 1 / (2 * (n_pop + 1))
	    max_val = 1 - min_val
	    af_pop[af_pop < min_val] = min_val
	    af_pop[af_pop > max_val] = max_val
	    af[:,i] = af_pop
	    del L_pop, af_pop
	  np.save(args.out + ".pop_af", af)
	  print("Saved reference population allele frequencies as " + str(args.out) + \
	       ".pop_af.npy (Binary - np.float32)\n")
	  print("Column order of populations is: " + str(pops))
	  
	  ##  Fisher information
	  print("Estimating Fisher information.")
	  f_obs, ne_obs = fisher.fisher_obs(L, af, IDs, args.threads)
	  np.save(args.out + ".fisher_obs", f_obs)
	  print("Saved reference population observed Fisher information as " + str(args.out) + \
	    ".fisher_obs.npy (Binary - np.float32)\n")
	  np.save(args.out + ".ne_obs", ne_obs)
	  print("Saved reference population effective sample size estimates as " + str(args.out) + \
	    ".ne_obs.npy (Binary - np.float32)\n")
	  print("Column order of populations is: " + str(pops))
	  
	  if args.ne_obs_ind:
	    print("Estimating individual effective sample sizes.")
	    ne_ind_full = fisher.fisher_obs_ind(L, af, IDs, args.threads)
	    np.savetxt(args.out + ".ne_ind_full.txt", ne_ind_full.reshape(-1,1), fmt="%.7f")
	    # ne_ind_df = np.hstack((IDs, ne_ind_full.reshape(-1,1)))
	    # np.savetxt(args.out + ".ne_obs_ind.txt", ne_ind_df)
	    print("Save individual effective sample sizes as " + str(args.out) + \
	        ".ne_ind_full.txt")
	    
	  if args.loo:
	    print("Performing leave-one-out cross validation.")
	    logl_mat_loo = glassy.loo(L, af, IDs, args.threads, args.maf_iter, args.maf_tole)
	    np.savetxt(args.out + ".pop_like_LOO.txt", logl_mat_loo, fmt="%.7f")
	    print("Save leave-one-out cross validation log likelihoods as " + str(args.out) + \
	         ".pop_like_LOO.txt")
	    print("Column order of populations is: " + str(pops))
	  del af

	# Population assignment likelihood
	if args.get_pop_like:
	  print("Parsing population allele frequency file.")
	  assert os.path.isfile(args.pop_af_file), "Population allele frequency file does not exist!!"
	  A = np.load(args.pop_af_file)
	  print("Calculating likelihood of population assignment")
	  logl_mat = glassy.assignLL(L, A, args.threads)
	  np.savetxt(args.out + ".pop_like.txt", logl_mat, fmt="%.7f")
	  print("Saved population assignment log likelihoods as " + str(args.out) + \
	       ".pop_like.txt (text)")
	############################################################################
	# Z-score
	if args.get_reference_z_score:
	  # read in data
	  # Reference pop IDs
	  print("Parsing population ID file.")
	  assert os.path.isfile(args.pop_af_IDs), "Population ID file does not exist!!"
	  IDs = np.loadtxt(args.pop_af_IDs, delimiter = "\t", dtype = "str")
	  # Reference pop allele frequencies
	  print("Parsing population allele frequency file.")
	  assert os.path.isfile(args.pop_af_file), "Population allele frequency file does not exist!!"
	  A = np.load(args.pop_af_file)
	  # Reference pop allele depths
	  print("Parsing individual allele depths file.")
	  assert os.path.isfile(args.ind_ad_file), "Individual allele depths file does not exist!"
	  AD = np.load(args.ind_ad_file)
	  
	  # Unique reference pop names
	  pops = np.unique(IDs[:,1])
	  # number of individuals from beagle
	  n = L.shape[1] // 2
	  # Check number of individuals from beagle is same as reference file
	  assert (n == IDs.shape[0]), "Number of individuals in beagle and reference ID file do not match!"
	  if args.allele_count_threshold is not None:
	    assert (args.allele_count_threshold >= 0), "Allele count threshold needs to be greater than/equal to 0!"
	    allele_count_threshold = args.allele_count_threshold
	  else:
	    allele_count_threshold = 0
	  if args.ind_start is not None:
	    assert (args.ind_start > 0 and args.ind_start <= n), "Start individual index needs to be within range of number of individuals!"
	    ind_start = args.ind_start
	  else:
	    ind_start = 0
	  if args.ind_end is not None:
	    assert (args.ind_end > 0 and args.ind_end <= n), "End individual index needs to be within range of number of individuals!"
	    ind_end = args.ind_end
	  else:
	    ind_end = n
	  n_sub = ind_end - ind_start
	  z_out = np.empty((n_sub, 1), dtype = np.float32)
	  for i in range(ind_start, ind_end):
	    i_key = IDs[i,1]
	    k = np.argwhere(pops == i_key)[0][0]
	    AD_summary_dict, AD_array = zscore.AD_summary(L, AD, i, allele_count_threshold, args.single_read_threshold)
	    L_keep, loci_kept = zscore.get_L_keep(L, AD, AD_summary_dict, AD_array, i)
	    AD_factorial, AD_like, AD_index = zscore.get_factorials(AD_array, AD_summary_dict, 0.01)
	    # leave one out allele frequency calculation for reference
	    pop_index = np.argwhere(IDs[:,1] == i_pop)
	    pop_index_except = pop_index[pop_index != i]
	    L1 = pop_index_except * 2
	    L2 = L1 + 1
	    L_cat = np.concatenate((L1, L2))
	    L_cat_index = np.sort(L_cat, axis = 0).reshape(-1)
	    L_pop = np.ascontiguousarray(L[L_keep,L_cat_index])
	    af_pop = emMAF.emMAF(L_pop, maf_iter, maf_tole, t)
	    n_pop = L_pop.shape[1] // 2
	    min_val = 1 / (2 * (n_pop + 1))
	    max_val = 1 - min_val
	    af_pop[af_pop < min_val] = min_val
	    af_pop[af_pop > max_val] = max_val
	    del L_pop
	    # continue with new allele frequencies
	    W_l_obs, W_l_array = zscore.get_expected_W_l(L, L_keep, af_pop, AD, AD_array, AD_factorial, AD_like, AD_index, args.threads, i, k)
	    var_W_l_array = zscore.get_var_W_l(L, L_keep, af_pop, AD, AD_array, AD_factorial, AD_like, AD_index, W_l_array, args.threads, i, k)
	    z_mu = np.sum(W_l_array)
	    z_var = np.sum(var_W_l_array)
	    z_tmp = (W_l_obs - z_mu) / np.sqrt(z_var)
	    print("Finished individual " + str(i))
	    # print(AD_factorial)
	    # print(AD_like)
	    # print(AD_index)
	    print("z_mu: " + str(z_mu))
	    print("z_var: " + str(z_var))
	    print("z_obs: " + str(W_l_obs))
	    print("Loci used: " + str(loci_kept))
	    print("Z-score: " + str(z_tmp))
	    z_out[i-ind_start,0] = z_tmp
	  np.savetxt(args.out + ".reference_z_ind.txt", z_out, fmt="%.7f")
	  print("Saved " + str(n_sub) + " individual z-scores as " + str(args.out) + \
	       ".reference_z_ind.txt (text)")
	
	if args.get_assignment_z_score:
	  # read in data
	  # Reference pop IDs
	  print("Parsing population ID file.")
	  assert os.path.isfile(args.pop_af_IDs), "Population ID file does not exist!!"
	  IDs = np.loadtxt(args.pop_af_IDs, delimiter = "\t", dtype = "str")
	  # Reference pop allele frequencies
	  print("Parsing population allele frequency file.")
	  assert os.path.isfile(args.pop_af_file), "Population allele frequency file does not exist!!"
	  A = np.load(args.pop_af_file)
	  # Reference pop allele depths
	  print("Parsing individual allele depths file.")
	  assert os.path.isfile(args.ind_ad_file), "Individual allele depths file does not exist!"
	  AD = np.load(args.ind_ad_file)
	  
	  # Unique reference pop names
	  pops = np.unique(IDs[:,1])
	  # number of individuals from beagle
	  n = L.shape[1] // 2
	  # Check number of individuals from beagle is same as reference file
	  assert (n == IDs.shape[0]), "Number of individuals in beagle and reference ID file do not match!"
	  if args.allele_count_threshold is not None:
	    assert (args.allele_count_threshold >= 0), "Allele count threshold needs to be greater than/equal to 0!"
	    allele_count_threshold = args.allele_count_threshold
	  else:
	    allele_count_threshold = 0
	  if args.ind_start is not None:
	    assert (args.ind_start > 0 and args.ind_start <= n), "Start individual index needs to be within range of number of individuals!"
	    ind_start = args.ind_start
	  else:
	    ind_start = 0
	  if args.ind_end is not None:
	    assert (args.ind_end > 0 and args.ind_end <= n), "End individual index needs to be within range of number of individuals!"
	    ind_end = args.ind_end
	  else:
	    ind_end = n
	  n_sub = ind_end - ind_start
	  z_out = np.empty((n_sub, 1), dtype = np.float32)
	  for i in range(ind_start, ind_end):
	    pop_key = IDs[i,1]
	    k = np.argwhere(pops == pop_key)[0][0]
	    AD_summary_dict, AD_array = zscore.AD_summary(L, AD, i, allele_count_threshold, args.single_read_threshold)
	    L_keep, loci_kept = zscore.get_L_keep(L, AD, AD_summary_dict, AD_array, i)
	    AD_factorial, AD_like, AD_index = zscore.get_factorials(AD_array, AD_summary_dict, 0.01)
	    W_l_obs, W_l_array = zscore.get_expected_W_l(L, L_keep, A, AD, AD_array, AD_factorial, AD_like, AD_index, args.threads, i, k)
	    var_W_l_array = zscore.get_var_W_l(L, L_keep, A, AD, AD_array, AD_factorial, AD_like, AD_index, W_l_array, args.threads, i, k)
	    z_mu = np.sum(W_l_array)
	    z_var = np.sum(var_W_l_array)
	    z_tmp = (W_l_obs - z_mu) / np.sqrt(z_var)
	    print("Finished individual " + str(i))
	    print("z_mu: " + str(z_mu))
	    print("z_var: " + str(z_var))
	    print("z_obs: " + str(W_l_obs))
	    print("Loci used: " + str(loci_kept))
	    print("Z-score: " + str(z_tmp))
	    z_out[i-ind_start,0] = z_tmp
	  np.savetxt(args.out + ".z_ind.txt", z_out, fmt="%.7f")
	  print("Saved " + str(n_sub) + " individual z-scores as " + str(args.out) + \
	       ".z_ind.txt (text)")
	  
	############################################################################  
	# Mixture proportions
	if args.get_em_mix:
	  print("Parsing population assignment likelihood file.")
	  assert os.path.isfile(args.pop_like), "Population assignment log likelihood file does not exist!!"
	  assert os.path.isfile(args.pop_like_IDs), "ID file does not exist!!"
	  logl_mat = np.loadtxt(args.pop_like)
	  logl_mat_index = np.loadtxt(args.pop_like_IDs, delimiter = "\t", dtype = "str")
	  print("Calculating mixture proportions with EM")
	  em_out = mixture.em_mix(logl_mat, logl_mat_index, args.mixture_iter)
	  np.savetxt(args.out + ".em_mix.txt", em_out, fmt="%s")
	  print("Saved EM mixture proportions " + str(args.out) + \
	       ".em_mix.txt (text)")

	if args.get_mcmc_mix:
	  print("Parsing population assignment likelihood file.")
	  assert os.path.isfile(args.pop_like), "Population assignment log likelihood file does not exist!!"
	  assert os.path.isfile(args.pop_like_IDs), "ID file does not exist!!"
	  logl_mat = np.loadtxt(args.pop_like)
	  logl_mat_index = np.loadtxt(args.pop_like_IDs, delimiter = "\t", dtype = "str")
	  print("Calculating mixture proportions with EM")
	  mcmc_out = mixture.mcmc_mix(logl_mat, logl_mat_index, args.mixture_iter)
	  np.savetxt(args.out + ".em_mix.txt", mcmc_out, fmt="%s")
	  print("Saved MCMC mixture proportions " + str(args.out) + \
	       ".mcmc_mix.txt (text)")
############################################################
##### Define main #####
if __name__ == "__main__":
	main()
