"""
Get allele counts of major and minor
"""

__author__ = "Matt DeSaix"

import numpy as np
import sys

raw_counts_file = sys.argv[1]
raw_counts = np.loadtxt(raw_counts_file, dtype="int", skiprows=1)

source_majmin_file = sys.argv[2]
source_majmin = np.loadtxt(source_majmin_file, dtype="int", skiprows=1, usecols=(1,2))

n_sites = raw_counts.shape[0]
n_ind = raw_counts.shape[1] // 4

majmin_counts = np.empty((n_sites, n_ind*2), dtype=np.int32)
count_ind_index = np.tile(np.repeat(np.arange(n_ind), 2), (n_sites, 1)) * 4
count_majmin_index = np.tile(source_majmin, n_ind)
count_index = count_ind_index + count_majmin_index
majmin_counts = np.take_along_axis(raw_counts, count_index, 1).astype(np.int32)

np.savetxt(raw_counts_file + ".majmin.counts.txt.gz", majmin_counts, fmt = "%d")
