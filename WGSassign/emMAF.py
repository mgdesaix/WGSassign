"""
EM MAF
"""

__author__ = "Matt DeSaix"

# Libraries
import numpy as np

# Import scripts
from WGSassign import emMAF_cy

##### Functions #####
### Estimate MAF ###
def emMAF(L, iter, tole, t):
    m = L.shape[0] # Number of sites
    f = np.empty(m, dtype=np.float32)
    f.fill(0.25) # Uniform initialization
    f_prev = np.copy(f)
    for i in range(iter):
        emMAF_cy.emMAF_update(L, f, t)
        diff = emMAF_cy.rmse1d(f, f_prev)
        if diff < tole:
            print("EM (MAF) converged at iteration: " + str(i+1))
            break
        f_prev = np.copy(f)
    return f
