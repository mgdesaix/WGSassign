# WGSassign
Population assignment methods for whole-genome sequence and genotype likelihood data

## Get WGSassign and build
### Dependencies
The WGSassign software relies on the following Python packages that you can install through [mamba](https://github.com/mamba-org/mamba)/[conda](https://docs.conda.io/projects/conda/en/latest/index.html) (recommended) or pip:

- numpy
- cython
- scipy

You can create an environment through conda easily as follows:
```
# Conda environment
conda env create -f environment.yml
```

## Install and build
```bash
git clone https://github.com/mgdesaix/WGSassign.git
cd WGSassign
python setup.py build_ext --inplace
pip3 install -e .
```

You can now run WGSassign with the `WGSassign` command.

**Warning:** Installation on Mac OS X may require additional steps to enable OpenMP support. We recommend installation on HPC systems to take advantage of the parallelized code.

## Usage
### Running WGSassign
WGSassign works directly on genotype likelihood files. Beagle genotype likelihood files can be generated from BAM files using [ANGSD](https://github.com/ANGSD/angsd). 

```bash
# Estimate reference population allele frequencies
# Output = reference.popAF.npy (numpy binary matrix of size L (# loci) rows x K (ref pops) columns)
WGSassign -b reference.beagle.gz --pop_af_IDs reference_IDS.txt --get_reference_af --o reference

# Estimate population assignment likelihoods
# Output = assign.pop_like.txt (text file of size N (individuals) rows x K (ref pops) columns)
WGSassign -b assign.beagle.gz --pop_af_file reference.popAF.npy --get_pop_like --o assign
```

## Acknowledgements

This software makes extensive use of python/cython functions from [PCAngsd](https://github.com/Rosemeis/pcangsd). Most notably, we have incorporated their parallelized cython functions for efficiently reading large Beagle files and also expectation-maximization estimation of allele frequencies. More details regarding the uses of PCangsd can be found in the following publications:

Population structure: [Inferring Population Structure and Admixture Proportions in Low-Depth NGS Data](http://www.genetics.org/content/210/2/719).\
HWE test: [Testing for Hardy‐Weinberg Equilibrium in Structured Populations using Genotype or Low‐Depth NGS Data](https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.13019).\
Selection: [Detecting Selection in Low-Coverage High-Throughput Sequencing Data using Principal Component Analysis](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04375-2).
