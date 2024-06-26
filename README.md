# RaMeDiES: Rare Mendelian Disorder Enrichment Statistics

This software package implements three well-calibrated statistical methods for analyzing cohorts of rare disease patients to find:
1. genes recurrently impacted by _de novo_ mutations across the cohort
2. genes recurrently impacted by inherited compound heterozygous variants across the cohort
3. genes harboring significant compound heterozygous variants in individual patients

Our [RaMeDiES wiki](https://github.com/hms-dbmi/RaMeDiES/wiki) also details how we ran our [pathway analysis](https://github.com/hms-dbmi/RaMeDiES/wiki/Pathway-analysis) to find pathways enriched with
candidate diagnostic variants across phenotypically similar patients.

If you use RaMeDiES in your work, please cite our publication: 
> SN Kobren*, MA Moldovan*, R Reimers, D Traviglia, X Li, D Barnum, A Veit, RI Corona, GdV Carvalho Neto, J Willett, M Berselli, W Ronchetti, SF Nelson, JA Martinez-Agosto, R Sherwood, J Krier, IS Kohane, Undiagnosed Diseases Network, SR Sunyaev (2024). "Joint, multifaceted genomic analysis enables diagnosis of diverse, ultra-rare monogenic presentations." _bioRxiv._ doi: [10.1101/2024.02.13.580158](https://www.biorxiv.org/content/10.1101/2024.02.13.580158).

## :sparkles: Prerequisites
* Python 3.6+, R 4.1+
* Python libraries: os, sys, argparse v1.1+, numpy v1.23.3+, scipy v1.91+, rpy2 v3.15.16+, requests v3.31+, urllib3 v1.26.8+
* R packages: cluster
* :exclamation: **Operating System:** Linux or MacOS;  Windows is not supported.
  
## :sparkles: Configuration
Edit the configuration `cfg.py` file to include the full path to your local installation of this repository.

```
script_directory = "/full/path/to/github/repo/RaMeDiES/"
```

## :sparkles: Precomputed data files
We have precomputed per-gene mutational targets for various variant functionality scores with respect to GRCh38/hg38. *The most up-to-date versions of these files can be found in* `/full/path/to/github/repo/RaMeDiES/data`.

A freeze of the precomputed files used in our initial manuscript submission (2024-02-01) can be downloaded from [Harvard Dataverse](https://doi.org/10.7910/DVN/UISZTE).

* `ens2gene.txt.gz` (136 KB)
* `pseudogenes.txt.gz` (231 KB)
* `score_lists_CI.txt.gz` (629 KB)
* `score_lists_CS.txt.gz` (97 MB)
* `score_lists_II.txt.gz` (218 KB)
* `score_lists_IS.txt.gz` (13.35 MB)
* `shet_table.txt.gz` (489 KB)

## :sparkles: RaMeDiES Framework
Descriptions of, sample code for running, and customizable parameters for the following steps of our statistical framework are detailed in [our wiki](https://github.com/hms-dbmi/RaMeDiES/wiki): 

* [Preprocess input variant data](https://github.com/hms-dbmi/RaMeDiES/wiki/Variant-data-input)
* [Cohort-level _de novo_ recurrence](https://github.com/hms-dbmi/RaMeDiES/wiki/Cohort-level-de-novo-recurrence)
* [Cohort-level compound heterozygosity](https://github.com/hms-dbmi/RaMeDiES/wiki/Cohort-level-compound-heterozygosity)
* [Individual-level compound heterozygosity](https://github.com/hms-dbmi/RaMeDiES/wiki/Individual-level-compound-heterozygosity)
* [Gene set enrichment across patient subgroups](https://github.com/hms-dbmi/RaMeDiES/wiki/Pathway-analysis)

:exclamation: Our statistical models operate at the level of "mutational targets" rather than individual-level variant data. These intermediate computed files can be shared freely to enable cross-cohort meta-analyses! See our [enabling cross-cohort analyses](https://github.com/hms-dbmi/RaMeDiES/wiki/Metaanalyses) wiki page for more details.

## :sparkles: Contact
If you have questions or comments about running any of the code found in this repository, please contact Shilpa Kobren or Mikhail Moldovan at [first name]_[last name] at hms.harvard.edu.
