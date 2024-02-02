# RaMeDiES: Rare Mendelian Disorder Enrichment Statistics

This software package implements three well-calibrated statistical methods for analyzing cohorts of rare disease patients to find:
1. genes recurrently impacted by _de novo_ mutations across the cohort
2. genes recurrently impacted by inherited compound heterozygous variants across the cohort
3. genes harboring significant compound heterozygous variants in individual patients

If you use RaMeDiES in your work, please cite our publication: 
> SN Kobren*, MA Moldovan*, R Reimers, D Traviglia, X Li, D Barnum, A Veit, J Willett, R Sherwood, J Krier, IS Kohane, Undiagnosed Diseases Network, SR Sunyaev (2024). "Joint, multifaceted genomic analysis enables diagnosis of diverse, ultra-rare monogenic presentations." _arXiv._ doi:[XXX](https://www.google.com/).

## :sparkles: Getting Started
### Prerequisites
* Python 3.6 (or above)
* Python libraries: os, sys, argparse, numpy, scipy

### Configuration
Edit the configuration `cfg.py` file to include the full path to your local installation of this repository. All scripts expect this variable to end with a forward slash `/`. 

```
script_directory = "/full/path/to/github/directory/ramedies/"
```

### Download precomputed data files
All RaMeDiES statistical models operate at the level of _mutational targets_, which intuitively correspond to the total mutation rate of all possible variants (of a particular type) within a gene (with as high a functionality score or higher). We have precomputed per-gene mutational targets for CADD and SpliceAI variant functionality scores with respect to GRCh38/hg38. You must download these files and store them locally in `/full/path/to/github/directory/ramedies//data`.

*Seven* required files must be downloaded from [Harvard Dataverse](https://doi.org/10.7910/DVN/UISZTE): 
* `ens2gene.txt.gz`
* `pseudogenes.txt.gz`
* `score_lists_CI.txt.gz`
* `score_lists_CS.txt.gz`
* `score_lists_II.txt.gz`
* `score_lists_IS.txt.gz`
* `shet_table.txt.gz`

## :sparkles: Format Input Variant Data
We expect annotated [VCF (variant call format)](https://samtools.github.io/hts-specs/VCFv4.2.pdf) input files with variants aligned to GRCh38/hg38 for use with RaMeDiES. The `INFO` field of the VCF file should include the following parameters, although the names of these columns can be specified separately in the `cfg.py` file. 

```
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  SAMPLE_GENOTYPE
20  14370  .  G  A  29  PASS  GeneID=, etc.  GT:DP  0/1:40
```
