"""
cfg.py

Configuration file for RaMeDiES algorithms
"""

import os

# Automatically extract current directory (where this script is)
script_directory = os.path.dirname(os.path.realpath(__file__))
script_directory = '/'.join(script_directory.split('/')[:-1])

# Paths to files containing mutation rate-score distributions
# Change the values in this dictionary to incorporate custom annotations.
# See the GitHub wiki for details.

variant_scores_files = {'C': {'S': {"CADD": f"{script_directory}/data/CADD_CS.txt.gz",
                                    "AlphaMissense": f"{script_directory}/data/AlphaMissense_MS.txt.gz",
                                    "REVEL": f"{script_directory}/data/REVEL_MS.txt.gz",
                                    "PAI3D": f"{script_directory}/data/PAI3D_MS.txt.gz"},
                              'I': f"{script_directory}/data/CADD_CI.txt.gz"},
                        'I': {'S': f"{script_directory}/data/SpliceAI_IS.txt.gz",
                              'I': f"{script_directory}/data/SpliceAI_II.txt.gz"}}

# Path to pseudogene/RNA gene/overlapping gene list
pseudogenes = f"{script_directory}/data/pseudogenes.txt.gz"

# Path to the table containing GeneBayes, s_het, and other gene constraint values
gene_score_table = f"{script_directory}/data/gene_constraint_scores.txt.gz"

# Path to the table with ENSEMBL ID-Gene ID matches
ens2gene = f"{script_directory}/data/ens2gene.txt.gz"

# Drop homozygous recessive variants bool
drop_HR = True

# Dictionary with VEP consequences
# 'U': UTR
# 'I': intronic
# 'S': synonymous
# 'C': coding
# 'O': other

# To incorporate custom annotation, change the keys of this dictionary while retaining the values.
# For the incorporation of other variant types (such as regulatory region variants), you will need to precompute the
#  variant score files and include them in the "variant_score_files" dictionary defined above. 
VEP_cons_dict = {'5_prime_UTR_variant': 'U',
                 '5_prime_UTR': 'U',
                 '5PRIME_UTR': 'U',
                 '3_prime_UTR_variant': 'U',
                 '3_prime_UTR': 'U',
                 '3PRIME_UTR': 'U',
                 'upstream_gene_variant': 'U',
                 'upstream': 'U',
                 'downstream_gene_variant': 'U',
                 'downstream': 'U',
                 'DOWNSTREAM': 'U',
                 'regulatory': 'U',
                 'UPSTREAM': 'U',
                 'intron_variant': 'I',
                 'intron': 'I',
                 'INTRONIC': 'I',
                 'splice_acceptor_variant': 'I',
                 'splice_acceptor': 'I',
                 'splice_donor_variant': 'I',
                 'splice_donor': 'I',
                 'splice_donor_region_variant': 'I',
                 'splice_region_variant': 'I',
                 'splice_donor_5th_base_variant': 'I',
                 'splice_donor_5th_base': 'I',
                 'splice_polypyrimidine_tract_variant': 'I',
                 'splice_polypyrimidine_tract': 'I',
                 'splice': 'I',
                 'CANONICAL_SPLICE': 'I',
                 'SPLICE_SITE': 'I',
                 'synonymous_variant': 'S',
                 'synonymous': 'S',
                 'stop_retained_variant': 'C',
                 'stop_lost': 'C',
                 'stop_gained': 'C',
                 'STOP_GAINED': 'C',
                 'start_lost': 'M',
                 'STOP_LOST': 'C',
                 'start_retained_variant': 'M',
                 'missense_variant': 'M',
                 'missense': 'M',
                 'inframe_deletion': 'M',
                 'inframe_insertion': 'M',
                 'frameshift_variant': 'M',
                 'FRAME_SHIFT': 'M',
                 'INFRAME': 'M',
                 'protein_altering_variant': 'M',
                 'incomplete_terminal_codon_variant': 'C',
                 'coding_sequence_variant': 'C',
                 'coding_transcript_variant': 'C',
                 'NON_SYNONYMOUS': 'C',
                 'NMD': 'C',
                 'regulatory_region': 'O',
                 'transcript_ablation': 'O',
                 'transcript_amplification': 'O',
                 'feature_elongation': 'O',
                 'feature_truncation': 'O',
                 'mature_miRNA_variant': 'O',
                 'non_coding_transcript_variant': 'O',
                 'non_coding_exon': 'O',
                 'non_coding': 'O',
                 'TFBS_ablation': 'O',
                 'TFBS_amplification': 'O',
                 'TF_binding_site_variant': 'O',
                 'regulatory_region_ablation': 'O',
                 'regulatory_region_amplification': 'O',
                 'regulatory_region_variant': 'O',
                 'intergenic_variant': 'O',
                 'intergenic': 'O',
                 'sequence_variant': 'O',
                 'NONCODING_CHANGE': 'O',
                 '': None}

# Dictionary specifying inheritance pattern keywords 
# Change the keys of this dictionary in case of other annotation scheme.
inherited_from_dict = {"mom": 'M',
                       "dad": 'P',
                       "": "DN",
                       "neither": "DN"}

# Dictionary specifying variant file format headers
# Change the values of this dictionary to incorporate custom formats
vcf_format_dict = {"chrom": "chromosome",
                   "position": "1-indexed_location",
                   "ref_al": "ref_allele",  # Reference allele
                   "alt_al": "alt_allele",  # Alternative allele
                   "var_annot": "consequence",  # Variant annotation
                   "SAI_AG": "SpliceAI_acceptor-gain-score",
                   "SAI_AL": "SpliceAI_acceptor-loss-score",
                   "SAI_DG": "SpliceAI_donor-gain-score",
                   "SAI_DL": "SpliceAI_donor-loss-score",
                   "coding_score": "CADD-raw",  # CADD score of a mutation
                   "ensembl_gene_id": "ensembl_gene_id",  # ENSEMBL gene ID
                   "MAF": "MAF",
                   "inherited_from": "inherited_from",  # Parent annotation
                   "qual_track": "DenovoMutationRate"}

# Reverse dictionary with input variant file format headers
rev_vcf_format_dict = {v: k for k, v in vcf_format_dict.items()}

# Quality track keywords
# Change the keys of this dictionary to incorporate custom formats
qual_value_dict = {"high": True,
                   "low": False}

# List of variant annotations
# 1st letter: 'C' for coding, 'I' for intronic
# 2nd letter: 'S' for SNP, 'I' for indel
var_annot_list = [('C', 'S'), ('C', 'I'), ('I', 'S'), ('I', 'I')]

# Identifiers used in thresholding of scores of input variants
#   see parse_variant_lib/get_score_val function
var_type_to_var_str = {('C', 'S'): 'C',
                       ('C', 'I'): 'CInd',
                       ('I', 'S'): 'I',
                       ('I', 'I'): 'I'}

# output file suffix: total number of DENOVO variants observed across cohort for each of 4 variant types
varcount_sums_DN = "denovo_variant_counts"

# output file suffix: total number of COMPHET variants observed across cohort for each of 16 variant types
mutnum_prod_CH = "comphet_variant_counts"

# output file suffix: per-gene DENOVO mutational targets (comma-separated per patient) 
muttargs_list_DN_ID = "denovo_mutational_targets"

# output file suffix: per-gene COMPHET mutational targets (comma-separated per patient) 
muttargs_CH_ID = "comphet_mutational_targets"

# output file suffix: per-patient COMPHET variant count distribution (for identifying outliers for QC purposes) 
mutnum_prod_dist_CH = "comphet_variant_product_distribution"

# output file suffix: per-patient variant count distribution (for identifying outliers for QC purposes) by inheritance
varcount_mask = "variant_distribution"  # number of variants by type -> number of patients with that count

# output file suffix: de novo recurrence across cohort
DN_result = "denovo_cohort_recurrence"

# output file suffix: compound heterozygous recurrence across cohort
CH_result = "comphet_cohort_recurrence"

# output file suffix: individual-level compound heterozygous results
CH_IND_result = "comphet_individual_level"

# Precision for the p-value calculated with an infinite sum
pval_precision = 1e-6

# Binomial probability of the false diagnosis rate being higher than the given value
false_diag_rate = 0.05

# Maximal number of values calculated in the 'infinite' sum
maxIHval = 1000

# Threshold for the Irwin-Hall parameter above which the Irwin-Hall distribution is approximated with Normal
IH_norm_approx_thr = 10

# Written on 01.23.2024 by Mikhail Moldovan, HMS DBMI
