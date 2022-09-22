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
parser.add_argument("-p", "--plink", metavar="FILE-PREFIX",
	help="Prefix PLINK files (.bed, .bim, .fam)")
parser.add_argument("-t", "--threads", metavar="INT", type=int, default=1,
	help="Number of threads")
parser.add_argument("-o", "--out", metavar="OUTPUT", default="pcangsd",
	help="Prefix for output files")

# parser.add_argument("--maf_iter", metavar="INT", type=int, default=200,
# 	help="Maximum iterations for minor allele frequencies estimation - EM (200)")
# parser.add_argument("--maf_tole", metavar="FLOAT", type=float, default=1e-4,
# 	help="Tolerance for minor allele frequencies estimation update - EM (1e-4)")
parser.add_argument("--iter", metavar="INT", type=int, default=100,
	help="Maximum iterations for estimation of individual allele frequencies (100)")
parser.add_argument("--tole", metavar="FLOAT", type=float, default=1e-5,
	help="Tolerance for update in estimation of individual allele frequencies (1e-5)")

#############################################################3
# Step 1) Reference population allele frequencies
parser.add_argument("--pop_af_IDs", metavar="FILE",
	help="Filepath to IDs for reference population beagle")
parser.add_argument("--get_reference_af", action="store_true", 
  help="Estimate allele frequencies for reference populations")

# Step 2) Estimate likelihoods of assignment
parser.add_argument("--pop_af_file", metavar="FILE",
	help="Filepath to reference population allele frequencies")
parser.add_argument("--get_pop_like", action="store_true", 
  help="Estimate log likelihood of individual assignment to each reference population")
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
	assert (args.beagle is not None) or (args.plink is not None), \
			"Please provide input data (args.beagle or args.plink)!"

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
	from WGSassign import reader_cy
	from WGSassign import shared
	# from pcangsd import covariance
	# from pcangsd import selection
	# from pcangsd import inbreed
	# from pcangsd import admixture
	# from pcangsd import tree
	from WGSassign import glassy

	# Parse data
	if args.beagle is not None:
		print("Parsing Beagle file.")
		assert os.path.isfile(args.beagle), "Beagle file doesn't exist!"
		L = reader_cy.readBeagle(args.beagle)
		m = L.shape[0]
		n = L.shape[1]//2
	else:
		print("Parsing PLINK files.")
		# assert args.filter is None, "Please perform sample filtering in PLINK!"
		# assert args.filterSites is None, "Please perform site filtering in PLINK!"
		assert os.path.isfile(args.plink + ".bed"), "Bed file doesn't exist!"
		assert os.path.isfile(args.plink + ".bim"), "Bim file doesn't exist!"
		assert os.path.isfile(args.plink + ".fam"), "Fam file doesn't exist!"

		# Count number of sites and samples
		m = extract_length(args.plink + ".bim")
		n = extract_length(args.plink + ".fam")
		with open(args.plink + ".bed", "rb") as bed:
			G = np.fromfile(bed, dtype=np.uint8, offset=3)
		G_len = ceil(n/4)
		G = G.reshape(m, G_len)
		L = np.zeros((m, 2*n), dtype=np.float32)
		reader_cy.convertBed(L, G, G_len, args.plink_error, m, n, args.threads)
	print("Loaded " + str(m) + " sites and " + str(n) + " individuals.")
	m_old = L.shape[0] # For future reference
	
####################################
  # Step 1) Reference population allele frequencies
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
	  f = np.empty((m, npops), dtype=np.float32)
	  # For each reference population, estimate the allele frequencies from the beagle file
	  for i in range(npops):
	    # get indices of which rows in ID file correspond to the given reference pop
	    pop_index = np.argwhere(IDs[:,1] == pops[i])
	    # convert indices to relevant column indices of "L" file in pcangsd (Remember Beagle file is converted to 2 cols per individual)
	    L1 = pop_index*2
	    L2 = L1 + 1
	    L_cat = np.concatenate((L1, L2))
	    L_cat_index = np.sort(L_cat, axis = 0)
	    L_cat_index.shape = (len(L_cat_index))
	    L_pop = np.ascontiguousarray(L[:,L_cat_index])
	    f_pop = shared.emMAF(L_pop, args.maf_iter, args.maf_tole, args.threads)
	    f[:,i] = f_pop
	    del L_pop, f_pop
	  np.save(args.out + ".popAF", f)
	  del f
	  print("Saved reference population allele frequencies as " + str(args.out) + \
	       ".popAF.npy (Binary - np.float32)\n")

	# Step 2) Population assignment likelihood
	if args.get_pop_like:
	  print("Parsing population allele frequency file.")
	  assert os.path.isfile(args.pop_af_file), "Population allele frequency file does not exist!!"
	  # Need to figure out cython reader for this
	  # A = reader_cy.readPopAF(args.pop_af_file)
	  A = np.load(args.pop_af_file)
	  print("Calculating likelihood of population assignment")
	  logl_mat = glassy.assignLL(L, A, args.threads)
	  np.savetxt(args.out + ".pop_like", logl_mat, fmt="%.7f")
	  print("Saved population assignment log likelihoods as " + str(args.out) + \
	       ".pop_like (text)")

############################################################
##### Define main #####
if __name__ == "__main__":
	main()
