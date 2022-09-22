# WGSassign
Population assignment methods for whole-genome sequence and genotype likelihood data


This software makes extensive use of python/cython functions from [PCAngsd](https://github.com/Rosemeis/pcangsd). Most notably, we have incorporated their parallelized cython functions for efficiently reading large Beagle files and also expectation-maximization estimation of allele frequencies. More details regarding the uses of PCangsd can be found in the following publications:

Population structure: [Inferring Population Structure and Admixture Proportions in Low-Depth NGS Data](http://www.genetics.org/content/210/2/719).\
HWE test: [Testing for Hardy‐Weinberg Equilibrium in Structured Populations using Genotype or Low‐Depth NGS Data](https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.13019).\
Selection: [Detecting Selection in Low-Coverage High-Throughput Sequencing Data using Principal Component Analysis](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04375-2).