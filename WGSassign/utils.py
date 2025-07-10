
__author__ = "Eric C. Anderson"

import numpy as np
import pandas as pd
import gzip


def print_sample_and_site_summary(sample_names, site_names):
    def preview(name_list):
        n = len(name_list)
        if n <= 4:
            return ", ".join(name_list)
        else:
            return ", ".join(name_list[:2]) + ", ..., " + ", ".join(name_list[-2:])

    print(f"sample_names: {len(sample_names)} samples total: {preview(sample_names)}")
    print(f"site_names: {len(site_names)} sites total: {preview(site_names)}")



def filter_sites_to_common(L, site_names, site_names_target):
    """
    Keep only rows from L and site_names where site_names are present in site_names_target.
    Assumes site_names and site_names_target are lists of strings.
    
    Returns:
        L_filtered, site_names_filtered
    """
    site_names = np.array(site_names)
    site_names_target = set(site_names_target)  # convert once for fast lookup

    mask = np.isin(site_names, list(site_names_target))
    num_filtered = np.sum(~mask)

    if num_filtered > 0:
        print(f"\tFiltered out {num_filtered} sites not present in the target site list.")

    L_filtered = L[mask, :]
    site_names_filtered = site_names[mask].tolist()

    return L_filtered, site_names_filtered






def write_ass_mats(
    filename,
    loglike_mat,
    sample_names,
    pop_names,
    partition_count=1,
    print_part_column=True,
    sample_locations=None,
    doing_LOO=False
):
    """
    Write assignment likelihood matrix to a tab-delimited file, optionally gzipped.

    Parameters:
        filename : str
            Output file name. If it ends in '.gz', output will be compressed.
        loglike_mat : np.ndarray
            Assignment log-likelihoods, shape (n_ind * partition_count, K)
        sample_names : list of str
            Sample names, length n_ind
        pop_names : list of str
            Population names (column headers), length K
        partition_count : int
            Number of site partitions (1 = whole genome).
        print_part_column : bool
            Whether to include 'data_part' column in output.
        sample_locations : list of str or None
            If provided, location or population for each sample.
        doing_LOO : bool
            If True, check that sample_locations are a subset of pop_names.
    """

    n_ind = len(sample_names)
    K = len(pop_names)
    expected_shape = (n_ind * partition_count, K)

    if loglike_mat.shape != expected_shape:
        raise ValueError(f"loglike_mat shape mismatch: expected {expected_shape}, got {loglike_mat.shape}")

    if not print_part_column and partition_count != 1:
        raise ValueError("print_part_column=False is only allowed if partition_count == 1")

    if sample_locations is not None:
        if len(sample_locations) != n_ind:
            raise ValueError("Length of sample_locations does not match sample_names")
        if doing_LOO and not set(sample_locations).issubset(set(pop_names)):
            raise ValueError("sample_locations contains values not in pop_names (required for LOO mode)")

    sample_col = np.repeat(sample_names, partition_count)
    data = {"sample": sample_col}

    if sample_locations is not None:
        location_col = np.repeat(sample_locations, partition_count)
        if doing_LOO is True:
            location_col_name = "source_pop"
        else:
            location_col_name = "location"
        data[location_col_name] = location_col


    if print_part_column:
        part_col = np.tile(np.arange(partition_count), n_ind)
        data["data_part"] = part_col

    df_ll = pd.DataFrame(loglike_mat, columns=pop_names)
    df = pd.concat([pd.DataFrame(data), df_ll], axis=1)

    # Handle gzip output if filename ends with .gz
    if filename.endswith(".gz"):
        with gzip.open(filename, "wt") as f:
            df.to_csv(f, sep="\t", index=False, float_format="%.6f")
    else:
        df.to_csv(filename, sep="\t", index=False, float_format="%.6f")

    print(f"Wrote assignment matrix to {filename}")





def partition_loglikes(per_site_ll, partition_count):
    """
    Partition and sum 1D per-site log-likelihoods into a 1D vector of partition sums.

    Parameters:
        per_site_ll : np.ndarray, shape (L_sites,)
            Per-site log-likelihoods of one individual at L_sites sites for one population.
        partition_count : int
            Number of site partitions.

    Returns:
        logl_parts : np.ndarray, shape (partition_count,)
            Summed log-likelihoods for each partition.
    """
    if per_site_ll.ndim != 1:
        raise ValueError("per_site_ll must be a 1D array")

    L_sites = per_site_ll.shape[0]
    partition_labels = np.arange(L_sites) % partition_count  # break into groups by modulus
    logl_parts = np.zeros(partition_count, dtype=np.float32)
    np.add.at(logl_parts, partition_labels, per_site_ll)

    return logl_parts

