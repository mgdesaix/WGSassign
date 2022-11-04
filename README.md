# WGSassign
Population assignment methods for whole-genome sequence and genotype likelihood data

WGSassign was developed to conduct standard population assignment methods using genotype likelihood data in Beagle format. This software facilitates the estimation of allele frequencies across reference populations, allows users to perform leave-one-out cross validation with reference populations, and assign individuals of unknown origin to reference populations. See below for descriptions of the functions, and check out my [amre-adaptation](https://github.com/mgdesaix/amre-adaptation) GitHub repository for broader descriptions of WGSassigned being used for studying migratory bird populations.

## Get WGSassign and build
### Dependencies
The WGSassign software relies on the following Python packages that you can install through [mamba](https://github.com/mamba-org/mamba)/[conda](https://docs.conda.io/projects/conda/en/latest/index.html) (recommended) or pip:

- numpy (<=1.22.3, current WGSassign functionality fails with >=1.22.4)
- cython
- scipy

You can create an environment through mamba easily as follows using the `environment.yml` file provided:
```
# Conda environment
mamba env create -f environment.yml
```

## Install and build
```bash
git clone https://github.com/mgdesaix/WGSassign.git
cd WGSassign
mamba env create -f environment.yml
conda activate WGSassign
python setup.py build_ext --inplace
pip3 install -e .
```

You can now run WGSassign with the `WGSassign` command.

**Note:** WGSassign has only been tested on Linux systems

## Usage
### Running WGSassign
WGSassign works directly on genotype likelihood files. Beagle genotype likelihood files can be generated from BAM files using [ANGSD](https://github.com/ANGSD/angsd). 

Options:

```bash
--% WGSassign
usage: WGSassign [-h] [-b FILE] [-p FILE-PREFIX] [-t INT] [-o OUTPUT] [--maf_iter INT] [--maf_tole FLOAT]
                 [--pop_af_IDs FILE] [--get_reference_af] [--pop_af_file FILE] [--get_pop_like] [--loo]
                 [--pop_like FILE] [--pop_like_IDs FILE] [--get_em_mix] [--get_mcmc_mix] [--mixture_iter INT]

options:
  -h, --help            show this help message and exit
  -b FILE, --beagle FILE
                        Filepath to genotype likelihoods in gzipped Beagle format from ANGSD
  -t INT, --threads INT
                        Number of threads
  -o OUTPUT, --out OUTPUT
                        Prefix for output files
  --maf_iter INT        Maximum iterations for minor allele frequencies estimation - EM (200)
  --maf_tole FLOAT      Tolerance for minor allele frequencies estimation update - EM (1e-4)
  --pop_af_IDs FILE     Filepath to IDs for reference population beagle
  --get_reference_af    Estimate allele frequencies for reference populations
  --pop_af_file FILE    Filepath to reference population allele frequencies
  --get_pop_like        Estimate log likelihood of individual assignment to each reference population
  --loo                 Perform leave-one-out cross validation
  --pop_like FILE       Filepath to population assignment log likelihood file
  --pop_like_IDs FILE   Filepath to IDs for population assignment log likelihood file
  --get_em_mix          Estimate mixture proportions with EM algorithm
  --get_mcmc_mix        Estimate mixture proportions with MCMC algorithm
  --mixture_iter INT    Maximum iterations mixture estimation - EM (200)
```

### Reference population allele frequencies

To create a file of allele frequencies across all reference populations is straightforward: provide 1) a beagle file and 2) a reference ID file. The reference ID file needs to be tab delimited and have 2 columns, with each row being an individual and the 2nd column specifying the reference population the individual belongs to. It is **essential** that the sample order in the ID file (row-wise) reflects the sample order in the Beagle file (column-wise); though header sample names in the beagle file don't have to match the first column of the ID file. I do this by just creating a reference file from whatever BAM list I used to create the beagle file.

The following code examples can be run using the data provided with WGSassign in the `data/` directory.

Specifying `--get_reference_af` produces the reference population allele file. It is a numpy binary with L (# loci) rows x K (reference population) columns:

```bash
data_dir=/home/mgdesaix/projects/WGSassign-practice/data
breeding_beagle=amre.breeding.ind85.ds_2x.sites-filter.top_50_each.beagle.gz
breeding_IDs=amre.breeding.ind85.reference_k5.IDs.txt
outname=amre.breeding.ind85.ds_2x.sites-filter.top_50_each
# Estimate reference population allele frequencies using 20 threads
# Output = ${outname}.popAF.npy (numpy binary matrix of size L (# loci) rows x K (ref pops) columns)
WGSassign --beagle ${data_dir}/${breeding_beagle} --pop_af_IDs ${data_dir}/${breeding_IDs} --get_reference_af --out ./out/breeding/${outname} --threads 20
```

**Note:** numpy binaries are produced without population headers so double check the WGSassign output which states the population order of the columns!

### Leave-one-out cross validation

Cross-validation is an important technique for determining assignment accuracy, but recalculating allele frequencies each time you remove an individual can be slow. Fortunately, WGSassign is fast enough to provide leave-one-out assignment in a reasonable amount of time even for large of beagle files (ex. 
~ 5 million SNPs and 180 individuals = 30 min; 600k SNPs and 80 individuals = <1 min). When producing reference population allele frequencies (as described above), you can add `--loo` to specify the calculation of leave-one-out assignment likelihoods for each individual to each reference population. For LOO assignment, the allele frequency of the reference population the individual belongs to is recalculated without the individual, and then the individual is assigned to the different reference populations. The log-likelihoods of LOO assignment are output in a text file with N (individual) rows x K (reference population) columns. 

```bash
# Get likelihoods for leave-one-out assignment within known reference populations using 20 threads
# Output = 1) ${outname}.popAF.npy, 2) ${outname}.pop_like_LOO.txt
WGSassign --beagle ${data_dir}/${breeding_beagle} --pop_af_IDs ${data_dir}/${breeding_IDs} --get_reference_af --loo --out ./out/breeding/${outname} --threads 20
```

### Assigning individuals of unknown origin

Specifying `--get_pop_like` will produce the population log-likelihoods of assignment for each individual in a beagle file (typically unknown origin) to a reference population allele frequency file in numpy binary format.

```bash
nonbreeding_beagle=amre.nonbreeding.ind34.ds_2x.sites-filter.top_50_each.beagle.gz
nonbreeding_IDs=amre.nonbreeding.ind34.site.IDs.txt
outname2=amre.nonbreeding.ind34.ds_2x.sites-filter.top_50_each
# Estimate population assignment likelihoods
# Output = ${outname2}.pop_like.txt (text file of size N (individuals) rows x K (ref pops) columns)
WGSassign --beagle ${data_dir}/${nonbreeding_beagle} --pop_af_file ./out/breeding/${outname}.popAF.npy --get_pop_like --out ./out/nonbreeding/${outname2} --threads 20
```

### Mixture proportions

The mixture proportions of the population of assigned individuals can also be estimated. This only requires two files: 1) the likelihoods of assignment of individuals to the source populations, and 2) a tab-delimited ID file of the individuals, in the same order as the individual likelihood file, with the second column designating the *population* to consider the individuals a part of (ex. the nonbreeding site birds were sampled from). 

```sh
# Estimate mixture with EM algorithm
# Output = 1) ${outname2}.em_mix.txt (text file with column 1 as the mixed populations, and the remaining columns giving the mixture proportions of the source populations)
WGSassign --pop_like ./out/nonbreeding/${outname2}.pop_like.txt --pop_like_IDs ${data_dir}/${nonbreeding_IDs} --get_em_mix --out ./out/nonbreeding/${outname2}
```

## Acknowledgements

This software makes extensive use of python/cython functions from [PCAngsd](https://github.com/Rosemeis/pcangsd). Most notably, we have incorporated their parallelized cython functions for efficiently reading large Beagle files and also expectation-maximization estimation of allele frequencies. More details regarding the uses of PCangsd can be found in the following publications:

Population structure: [Inferring Population Structure and Admixture Proportions in Low-Depth NGS Data](http://www.genetics.org/content/210/2/719).\
HWE test: [Testing for Hardy‐Weinberg Equilibrium in Structured Populations using Genotype or Low‐Depth NGS Data](https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.13019).\
Selection: [Detecting Selection in Low-Coverage High-Throughput Sequencing Data using Principal Component Analysis](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04375-2).
