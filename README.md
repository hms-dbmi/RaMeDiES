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
Edit the configuration `cfg.py` file to specify the name of the column in your [input variant files](https://github.com/hms-dbmi/RaMeDiES/wiki/Variant-data-input) corresponding to the deleteriousness score you'd like to use for coding SNVs:

```python
vcf_format_dict = {"coding_snv_score" : "CADD-raw"}  # default is CADD-raw
```

Note that CADD is the only currently supported score for coding indels, and SpliceAI is the only currently supported score for intronic SNVs and indels.

## :sparkles: Precomputed data files
We have precomputed per-gene mutational targets for various variant functionality scores with respect to GRCh38/hg38. *Pointers to the most up-to-date versions of these files can be found in* `/full/path/to/github/repo/RaMeDiES/data`. 

You will need to have the [git Large File Storage (lfs)](https://git-lfs.com/) extension installed to download these files. It is supported on Mac, Windows, and Linux. On a Linux server without root access, you can install using: 

```bash
wget -O git-lfs.tar.gz https://github.com/git-lfs/git-lfs/releases/download/v3.5.1/git-lfs-linux-amd64-v3.5.1.tar.gz
tar -xvzf git-lfs.tar.gz
mv git-lfs/git-lfs ~/bin
git lfs install --skip-repo
```

Then you can run `git pull` to download the following files, or, if that fails, `git reset --hard` instead:

* `ens2gene.txt.gz` (133 KB)
* `pseudogenes.txt.gz` (225 KB)
* `gene_constraint_scores.txt.gz` (245 KB)
* `CADD_CI.txt.gz` (614 KB)
* `CADD_CS.txt.gz` (93 MB)
* `SpliceAI_II.txt.gz` (213 KB)
* `SpliceAI_IS.txt.gz` (13 MB)
* `AlphaMissense_MS.txt.gz` (68 MB)
* `PAI3D_MS.txt.gz` (63 MB)
* `REVEL_MS.txt.gz` (67 MB)
  
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
