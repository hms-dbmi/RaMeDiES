# RaMeDiES: Rare Mendelian Disorder Enrichment Statistics

This software package implements three well-calibrated statistical methods for analyzing cohorts of rare disease patients to find:
1. genes recurrently impacted by _de novo_ mutations across the cohort
2. genes recurrently impacted by inherited compound heterozygous variants across the cohort
3. genes harboring significant compound heterozygous variants in individual patients

If you use RaMeDiES in your work, please cite our publication: 
> SN Kobren*, MA Moldovan*, R Reimers, D Traviglia, X Li, D Barnum, A Veit, J Willett, M Berselli, W Ronchetti, R Sherwood, J Krier, IS Kohane, Undiagnosed Diseases Network, SR Sunyaev (2024). "Joint, multifaceted genomic analysis enables diagnosis of diverse, ultra-rare monogenic presentations." _bioRxiv._ doi: [10.1101/2024.02.13.580158](https://www.biorxiv.org/content/10.1101/2024.02.13.580158v1).

## :sparkles: Prerequisites
* Python 3.6+
* Python libraries: os, sys, argparse v1.1+, numpy v1.23.3+, scipy v1.91+
* :exclamation: **Operating System:** Linux distribution; compatibility on MacOS is not guaranteed, and Windows is not supported.
  
## :sparkles: Configuration
Edit the configuration `cfg.py` file to include the full path to your local installation of this repository. All scripts expect this variable to end with a forward slash `/`. 

```
script_directory = "/full/path/to/github/repo/RaMeDiES/"
```

## :sparkles: Download precomputed data files
All RaMeDiES statistical models operate at the level of _mutational targets_, which intuitively correspond to the total mutation rate of all possible variants (of a particular type) within a gene. We have precomputed per-gene mutational targets for CADD and SpliceAI variant functionality scores with respect to GRCh38/hg38. 

You must download these **seven** required files from [Harvard Dataverse](https://doi.org/10.7910/DVN/UISZTE) and store them locally in `/full/path/to/github/repo/RaMeDiES/data`:

* `ens2gene.txt.gz`
* `pseudogenes.txt.gz`
* `score_lists_CI.txt.gz`
* `score_lists_CS.txt.gz`
* `score_lists_II.txt.gz`
* `score_lists_IS.txt.gz`
* `shet_table.txt.gz`

## :sparkles: RaMeDiES Framework
Descriptions of, sample code for running, and customizable parameters for the following steps of our statistical framework are detailed in [our wiki](https://github.com/hms-dbmi/RaMeDiES/wiki): 

* [Preprocess input variant data](https://github.com/hms-dbmi/RaMeDiES/wiki/Variant-data-input)
* [Cohort-level _de novo_ recurrence](https://github.com/hms-dbmi/RaMeDiES/wiki/Cohort-level-de-novo-recurrence)
* [Cohort-level compound heterozygosity](https://github.com/hms-dbmi/RaMeDiES/wiki/Cohort-level-compound-heterozygosity)
* [Individual-level compound heterozygosity](https://github.com/hms-dbmi/RaMeDiES/wiki/Individual-level-compound-heterozygosity)

:exclamation: Our statistical models operate at the level of "mutational targets" rather than individual-level variant data. These intermediate computed files can be shared freely to enable cross-cohort meta-analyses! See our [enabling cross-cohort analyses](https://github.com/hms-dbmi/RaMeDiES/wiki/Metaanalyses) wiki page for more details.

## :sparkles: Contact
If you have questions or comments about running any of the code found in this repository, please contact Shilpa Kobren or Mikhail Moldovan at [first name]_[last name] at hms.harvard.edu.
