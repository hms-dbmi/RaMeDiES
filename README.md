# RaMeDiES: Rare Mendelian Disorder Enrichment Statistics

This software package implements three well-calibrated statistical methods for analyzing cohorts of rare disease patients to find:
1. genes recurrently impacted by _de novo_ mutations across the cohort
2. genes recurrently impacted by inherited compound heterozygous variants across the cohort
3. genes harboring significant compound heterozygous variants in individual patients

If you use RaMeDiES in your work, please cite our publication: 
> SN Kobren*, MA Moldovan*, R Reimers, D Traviglia, X Li, D Barnum, A Veit, J Willett, R Sherwood, J Krier, IS Kohane, Undiagnosed Diseases Network, SR Sunyaev (2024). "Joint, multifaceted genomic analysis enables diagnosis of diverse, ultra-rare monogenic presentations." _arXiv._ doi:[XXX](https://www.google.com/).

## :sparkles: Prerequisites
* Python 3.6 (or above)
* Python libraries: os, sys, argparse, numpy, scipy

## :sparkles: Configuration
Edit the configuration `cfg.py` file to include the full path to your local installation of this repository. All scripts expect this variable to end with a forward slash `/`. 

```
script_directory = "/full/path/to/github/directory/ramedies/"
```

## :sparkles: Download precomputed data files
All RaMeDiES statistical models operate at the level of _mutational targets_, which intuitively correspond to the total mutation rate of all possible variants (of a particular type) within a gene (with as high a functionality score or higher). We have precomputed per-gene mutational targets for CADD and SpliceAI variant functionality scores with respect to GRCh38/hg38. You must download these files and store them locally in `/full/path/to/github/directory/ramedies/data`.

*Seven* required files must be downloaded from [Harvard Dataverse](https://doi.org/10.7910/DVN/UISZTE): 
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
   * [Custom configuration of input variant data](https://github.com/hms-dbmi/RaMeDiES/wiki/Custom-configuration)    
* [Cohort-level _de novo_ recurrence](https://github.com/hms-dbmi/RaMeDiES/wiki/Cohort-level-de-novo-recurrence)
* [Cohort-level compound heterozygosity](https://github.com/hms-dbmi/RaMeDiES/wiki/Cohort-level-compound-heterozygosity)
* [Individual-level compound heterozygosity](https://github.com/hms-dbmi/RaMeDiES/wiki/Individual-level-compound-heterozygosity)

:exclamation: Our statistical models operate at the level of "mutational targets" rather than individual-level variant data. These intermediate computed files can be shared freely to enable cross-cohort meta-analyses! See our [enabling cross-cohort analyses](https://github.com/hms-dbmi/RaMeDiES/wiki/Metaanalyses) wiki page for more details.

## :sparkles: Contact
If you have questions or comments about running any of the code found in this repository, please contact Mikhail Moldovan or Shilpa Kobren at [first name]_[last name] at hms.harvard.edu.
