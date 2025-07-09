
__author__ = "Eric C. Anderson"

import numpy as np


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
