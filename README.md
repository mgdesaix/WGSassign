# WGSassign
Population assignment methods for whole-genome sequence and genotype likelihood data

WGSassign was developed to conduct standard population assignment methods using genotype likelihood data in Beagle format. This software facilitates the estimation of allele frequencies across reference populations, allows users to perform leave-one-out cross validation with reference populations, and assign individuals of unknown origin to reference populations. See below for descriptions of the functions, and check out my [amre-adaptation](https://github.com/mgdesaix/amre-adaptation) GitHub repository for broader descriptions of WGSassign being used for studying migratory bird populations.

## Get WGSassign and build
### Dependencies
The WGSassign software relies on the following Python packages that you can install through [mamba](https://github.com/mamba-org/mamba)/[conda](https://docs.conda.io/projects/conda/en/latest/index.html) (recommended) or pip:

- numpy (<=1.22.3, current WGSassign functionality fails with >=1.22.4)
- cython
- scipy

You can create an environment through conda/mamba easily as follows using the `environment.yml` file provided within WGSassign (see below).

### Install and build

Installation of WGSassign is simply:

```bash
git clone https://github.com/mgdesaix/WGSassign.git
cd WGSassign
mamba create -n WGSassign numpy scipy cython
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

```
--% WGSassign
usage: WGSassign [-h] [-b FILE] [-t INT] [-o OUTPUT] [--maf_iter INT] [--maf_tole FLOAT] [--pop_af_IDs FILE]
                 [--get_reference_af] [--pop_names FILE] [--ne_obs] [--loo] [--pop_af_file FILE] [--get_pop_like]
                 [--get_assignment_z_score] [--get_reference_z_score] [--ind_ad_file FILE] [--allele_count_threshold INT]
                 [--single_read_threshold] [--ind_start INT] [--ind_end INT] [--pop_like FILE] [--pop_like_IDs FILE]
                 [--get_em_mix] [--get_mcmc_mix] [--mixture_iter INT]

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
  --pop_af_IDs FILE     Filepath to individual IDs and populations for beagle
  --get_reference_af    Estimate allele frequencies for reference populations
  --pop_names FILE      Filepath to population names of allele frequency file
  --ne_obs              Estimate population and individuals effective sample sizes
  --loo                 Perform leave-one-out cross validation
  --pop_af_file FILE    Filepath to reference population allele frequencies
  --get_pop_like        Estimate log likelihood of individual assignment to each reference population
  --get_assignment_z_score
                        Calculate z-score for individuals
  --get_reference_z_score
                        Calculate z-score for individuals
  --ind_ad_file FILE    Filepath to individual allele depths
  --allele_count_threshold INT
                        Minimum number of loci needed to keep a specific allele count combination
  --single_read_threshold
                        Use only loci with a single read. Helpful for computational efficiency when individuals's sequencing depths vary.
  --ind_start INT       Start analysis at this individual index (0-index: i.e. 0 starts with the 1st individual)
  --ind_end INT         End analysis at this individual index (0-index: i.e. If you have 10 individuals, 9 is the 10th
                        individual)
  --pop_like FILE       Filepath to population assignment log likelihood file
  --pop_like_IDs FILE   Filepath to IDs for population assignment log likelihood file
  --get_em_mix          Estimate mixture proportions with EM algorithm
  --get_mcmc_mix        Estimate mixture proportions with MCMC algorithm
  --mixture_iter INT    Maximum iterations mixture estimation - EM (200)
```

## Reference population allele frequencies

To create a file of allele frequencies across all reference populations is straightforward: provide 1) a beagle file and 2) a reference ID file. The reference ID file needs to be tab delimited and have 2 columns, with each row being an individual and the 2nd column specifying the reference population the individual belongs to. It is **essential** that the sample order in the ID file (row-wise) reflects the sample order in the Beagle file (column-wise); though header sample names in the beagle file don't have to match the first column of the ID file. The easiest way to keep the proper order is to create a reference file from whatever BAM list you used to create the beagle file in ANGSD.

Allele frequency is parallelized across loci, so feel free to crank up the threads with `--threads`.

Specifying `--get_reference_af` produces the reference population allele file, `.pop_af.npy`. It is a numpy binary with L (# loci) rows x K (reference population) columns.

The following code examples can be run using the data provided with WGSassign in the `data/` directory. Note that population assignment with bird data is typically such that "breeding" individuals are of known origin and "nonbreeding" (e.g. wintering) individuals are of unknown breeding origin. Now you know why the following example variables and files are named such.

```bash
breeding_beagle=./data/amre.breeding.ind85.ds_2x.sites-filter.top_50_each.beagle.gz
breeding_IDs=./data/amre.breeding.ind85.reference_k5.IDs.txt
outname=./out/amre.breeding.ind85.ds_2x.sites-filter.top_50_each
# Estimate reference population allele frequencies using 20 threads
# Output = ${outname}.pop_af.npy (numpy binary matrix of size L (# loci) rows x K (ref pops) columns)
WGSassign --beagle ${breeding_beagle} --pop_af_IDs ${breeding_IDs} --get_reference_af --out ${outname} --threads 20
```

**Note:** Numpy binaries are produced without text headers. WGSassign will output a text file `.pop_names.txt` that specifies the order of the populations provided in the `.pop_af.npy` file.

### Reference effective sample sizes

When producing the reference population allele frequency file, WGSassign also has other features you can specify. `--ne_obs` adds the calculation of the observed Fisher information across loci and the output is a numpy binary of `.fisher_obs.npy` of the same *L* x *K* dimensions. The effective sample size of the reference populations for each locus are also calculated and the output is a numpy binary `.ne_obs.npy` of the same dimensions. The mean effective sample sizes of the reference populations is provided as output in a text file, `ne_obs.txt`, with the first row being the population names and the second row being the effective sample size of the population. `.ne_ind.txt` is a single column of each individual's effective sample size (same order as the file you provide to `--pop_af_ID`). 

```bash
WGSassign --beagle ${breeding_beagle} --pop_af_IDs ${breeding_IDs} --get_reference_af --ne_obs --out ${outname} --threads 20
```

### Leave-one-out cross validation

Cross-validation is an important technique for determining assignment accuracy, but recalculating allele frequencies each time you remove an individual can be slow. Fortunately, WGSassign is fast enough to provide leave-one-out assignment in a reasonable amount of time even for large of beagle files (ex. 
~ 5 million SNPs and 180 individuals = 30 min; 600k SNPs and 80 individuals = <1 min). When producing reference population allele frequencies (as described above), you can add `--loo` to specify the calculation of leave-one-out assignment likelihoods for each individual to each reference population. For LOO assignment, the allele frequency of the reference population the individual belongs to is recalculated without the individual, and then the individual is assigned to the different reference populations. The log-likelihoods of LOO assignment are output in a text file with N (individual) rows x K (reference population) columns. 

```bash
WGSassign --beagle ${breeding_beagle} --pop_af_IDs ${breeding_IDs} --get_reference_af --loo --out ${outname} --threads 20
```

## Assigning individuals of unknown origin

Specifying `--get_pop_like` will produce the population log-likelihoods of assignment (i.e. the likelihood of an individual's genotype arising from the given population allele frequency) for each individual in a provided beagle file to a reference population allele frequency file (numpy binary file; see previous code examples). The output text file is of size N (individuals) rows x K (reference population) columns.

```bash
nonbreeding_beagle=./data/amre.nonbreeding.ind34.ds_2x.sites-filter.top_50_each.beagle.gz
breeding_af=./out/amre.breeding.ind85.ds_2x.sites-filter.top_50_each.pop_af.npy
outname=./out/amre.nonbreeding.ind34.ds_2x.sites-filter.top_50_each
WGSassign --beagle ${nonbreeding_beagle} --pop_af_file ${breeding_af} --get_pop_like --out ${outname} --threads 20
```

You can then pull the text file into R or python and add your individual meta data file (remember, the output of the order of individuals in the rows is the same as the column order of your individuals in your beagle file!). It's then straightforward to calculate posterior probabilities of assignment and determine the population of highest likelihood of origin (heads up: if you plan on doing a z-score calculation for these individuals. Go ahead and create a nice, 2 column tab-delimited reference file with the individual samples names in the first column and the assigned population in the second column for later use). 

### Z-score

In population assignment, an individual will always have a population with the highest likelihood of producing that individual's genotype. However, it may be the case that the individual is actually from an unsampled population. But don't worry, you can test for this with WGSassign's implementation of a z-score of assignment! The general gist of the z-score is that you compare the observed likelihood of the individual being from a given population with the *expected* likelihood that would arise from an individual with the same read depth actually being from a given population. Whewf!

To do this, you'll need to get the allele counts for the major and minor alleles for all loci that are used to calculate allele frequency. The output you need for WGSassign is a numpy binary (`.npy` file), with L (# of loci) x 2 * N (individuals) columns. Where each individual has the major and minor allele count in their 2 columns. Produce this however you want, but here I'll explain how you can use ANGSD to do do this in a couple steps. 

**Get allele counts from ANGSD**

First, when producing the beagle files, also specify `-doCounts 1 -dumpCounts 4` to output the allele count data for all individuals, all alleles, across all loci (Sadly they don't just have a way to get major and minor alleles). But, it's then pretty straightforward to reference the major and minor allele columns from the beagle file (which conveniently enough are integers) to reduce this to the L x 2*N size file we actually need. You can use the script provided with WGSassign to create the allele count file: 

```sh
counts=your_counts_file.counts.gz
beagle=your_beagle_file.beagle.gz
./allele_counts_beagle.py ${counts} ${beagle}
```

The output is a numpy binary, `.majmin.counts.npy`, that can then be used within WGSassign.

**Get reference individuals' z-scores**

Ideally, z-scores for assignment of individuals to populations would produce a mean of 0 with a standard deviation of 1, but there's lots of variation. So, we will use the reference individuals to create a "standardized z-score" that can help shift individuals's z-scores (with unknown origin and hence we're using a z-score) to that mean of 0. 

Using the `${breeding_beagle}` and `${breeding_IDs}` files from before, and also the `pop_names.txt` file from `--get_reference_af`, along with our newly calculated numpy binary of allele counts as input to `--ind_ad_file`, we can get the reference z-scores. `--get_reference_z_score` is what allows us to do this.

```sh
WGSassign --beagle ${breeding_beagle} --pop_af_IDs ${breeding_IDs} --pop_names ${pop_names} --ind_ad_file ${breeding_ad} --get_reference_z_score --out ${outname} --threads 20
```

This outputs a text file of a single column with N (# of individuals) rows for each of their z-scores to their known population of origin. 

**Note:** Since z-score calculation is computationally intensive, you can speed things up somewhat, depending on the type of low-coverage data you have, by specifying `--single_read_threshold` which will limit the analysis to loci with one read (i.e. 1,0 or 0,1).

**Get unknown origin individuals' z-scores**

Calculating the z-score for individuals of unknown origin, requires that you have already used the `--get_pop_like` function to assign individuals. As noted in that section, you will need a file for input in `--pop_af_IDs` that is 2 columns, tab-delimited, with column one being the individuals sample names (same order as the beagle file), and column 2 having the assigned population names. The assigned population names have to match what you used for the reference population allele frequencies, which are also being input under `--pop_af_file`, and you are inputting the population names that were already produced as well, `--pop_names`. The allele count data for these individuals is also input, `--ind_ad_file`. Finally, you specify `--get_assignment_z_score` to get the z-scores, and as mentioned previously, you are welcome to add `--single_read_threshold`.

```sh
WGSassign --beagle ${nonbreeding_beagle} --pop_af_IDs ${nonbreeding_IDs} --pop_af_file ${breeding_af} --pop_names ${pop_names} --ind_ad_file ${nonbreeding_ad} --get_assignment_z_score --single_read_threshold --out ${outname} --threads 12
```

You can then pull the output z-score data from this, along with the reference z-scores, into R or whatever, and calculate the mean and standard deviation of the reference z-scores, and then use those to standardize the z-scores from the individuals of unkown origin.

## Mixture proportions (STILL IN DEVELOPMENT)

The mixture proportions of the population of assigned individuals can also be estimated. This only requires two files: 1) the likelihoods of assignment of individuals to the source populations, and 2) a tab-delimited ID file of the individuals, in the same order as the individual likelihood file, with the second column designating the *population* to consider the individuals a part of (ex. the nonbreeding site birds were sampled from). 

```sh
# Estimate mixture with EM algorithm
# Output = 1) ${outname2}.em_mix.txt (text file with column 1 as the mixed populations, and the remaining columns giving the mixture proportions of the source populations)
WGSassign --pop_like ./out/nonbreeding/${outname2}.pop_like.txt --pop_like_IDs ${data_dir}/${nonbreeding_IDs} --get_em_mix --out ./out/nonbreeding/${outname2}
```

## Acknowledgements

This study was funded by a Cooperative Agreement with the Alaska Department of Fish and Game (23-011) and an NSF CAREER award (008933-00002) to KCR. This work utilized the Alpine high performance computing resource at the University of Colorado Boulder. Alpine is jointly funded by the University of Colorado Boulder, the University of Colorado Anschutz, Colorado State University, and the National Science Foundation (award 2201538).

This software makes use of the code from [PCAngsd](https://github.com/Rosemeis/pcangsd) for efficiently reading Beagle files and also expectation-maximization estimation of allele frequencies (thank you!). More details regarding the uses of PCangsd can be found in the following publications:

Population structure: [Inferring Population Structure and Admixture Proportions in Low-Depth NGS Data](http://www.genetics.org/content/210/2/719).\
HWE test: [Testing for Hardy‐Weinberg Equilibrium in Structured Populations using Genotype or Low‐Depth NGS Data](https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.13019).\
Selection: [Detecting Selection in Low-Coverage High-Throughput Sequencing Data using Principal Component Analysis](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04375-2).
