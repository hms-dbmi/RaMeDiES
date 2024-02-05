"""
cfg

config file for ramediesDN

Before using the program, specify the directory in which the tool is located.
Directory name must end with "/ramedies". Variable script_directory
"""

# Change this variable to the directory in which the program is located
# Directory name must end with "/ramedies"
script_directory = "/n/data2/hms/dbmi/sunyaev/udn/workspace/mmoldovan/ramedies/"

# Paths to files containing mutation rate-score distributions
# Change the values in this dictionary to incorporate custom annotations.
# See the manual for details.
variant_scores_files = {'C': {'S' : f"{script_directory}data/score_lists_CS.txt.gz",
							  'I' : f"{script_directory}data/score_lists_CI.txt.gz"},
						'I': {'S' : f"{script_directory}data/score_lists_IS.txt.gz",
							  'I' : f"{script_directory}data/score_lists_II.txt.gz"}}

# Path to pseudogene/RNA gene/overlapping gene list
pseudogenes = f"{script_directory}/data/pseudogenes.txt.gz"

# Path to the table containing s_het values
shet_table = f"{script_directory}/data/shet_table.txt.gz"

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
# To incorporate custom annotation, change the keys of this dictionary
# while retaining the values.
# For the incorporation of other variant types such as regulatory region
# variants, calculate variant score files and include them into the 
# variant_scores_files. See the manual for details.
VEP_cons_dict = {'5_prime_UTR_variant' : 'U',
				 '3_prime_UTR_variant' : 'U',
				 'upstream_gene_variant' : 'U',
				 'downstream_gene_variant' : 'U',
				 'intron_variant' : 'I', 
				 'splice_acceptor_variant' : 'I',
				 'splice_donor_variant' : 'I',
				 'splice_donor_region_variant' : 'I',
				 'splice_region_variant' : 'I',
				 'splice_donor_5th_base_variant' : 'I',
				 'splice_polypyrimidine_tract_variant' : 'I',
				 'synonymous_variant' : 'S',
				 'stop_retained_variant' : 'C', 
				 'stop_lost' : 'C',
				 'stop_gained' : 'C',
				 'start_lost' : 'C', 
				 'start_retained_variant' : 'C', 
				 'missense_variant' : 'C', 
				 'inframe_deletion' : 'C',
				 'inframe_insertion' : 'C', 
				 'frameshift_variant' : 'C',
				 'protein_altering_variant' : 'C',
				 'incomplete_terminal_codon_variant' : 'C',
				 'coding_sequence_variant' : 'C',
				 'coding_transcript_variant' : 'C',
				 'frameshift_variant' : 'C',
				 'inframe_deletion' : 'C',
				 'inframe_insertion' : 'C',
				 'regulatory_region' : 'O',
				 'transcript_ablation' : 'O',
				 'transcript_amplification' : 'O',
				 'feature_elongation' : 'O',
				 'feature_truncation' : 'O',
				 'mature_miRNA_variant' : 'O',
				 'non_coding_transcript_variant' : 'O',
				 'TFBS_ablation' : 'O',
				 'TFBS_amplification' : 'O',
				 'TF_binding_site_variant' : 'O',
				 'regulatory_region_ablation' : 'O',
				 'regulatory_region_amplification' : 'O',
				 'regulatory_region_variant' : 'O',
				 'intergenic_variant' : 'O',
				 'sequence_variant' : 'O',
				 '' : None}

# Dictionary specifying inheritance pattern keywords 
# M: maternal
# P: paternal
# DN: de novo
# Change the keys of this dictionary in case of other annotation scheme.
inherited_from_dict = {"mom" : 'M', 
					   "dad" : 'P', 
					   "" : "DN",
					   "neither" : "DN"}

# Dictionary specifying VCF format headers
# Change the values of this dictionary to incorporate custom formats
vcf_format_dict = {"chrom" : "chromosome", 
				   "position" : "1-indexed_location", 
				   "ref_al" : "ref_allele", # Reference allele
				   "alt_al" : "alt_allele", # Alternative allele
				   "var_annot" : "consequence", # Variant annotation
				   "SAI_AG" : "SpliceAI_acceptor-gain-score", 
				   "SAI_AL" : "SpliceAI_acceptor-loss-score",
				   "SAI_DG" : "SpliceAI_donor-gain-score",
				   "SAI_DL" : "SpliceAI_donor-loss-score",
				   "CADD" : "CADD-raw", # CADD score of a mutation
				   "ENS_ID" : "ensembl_gene_id", # ENSEMBL gene ID
				   "MAF" : "MAF",
				   "inherited_from" : "inherited_from", # Parent annotation
				   "qual_track" : "DenovoMutationRate"} 

# Reverse dictionary with VCF format headers
rev_vcf_format_dict = {v : k for k, v in vcf_format_dict.items()}

# Quality track keywords
# Change the keys of this dictionary to incorporate custom formats
qual_value_dict = {"high" : True,
				   "low" : False}

# List of variant annotations
# 1st letter: 'C' for coding, 'I' for intronic
# 2nd letter: 'S' for SNP, 'I' for indel
var_annot_list = [('C' , 'S'), ('C', 'I'), ('I' , 'S'), ('I', 'I')]

# Identifier for the file with variant numbers for each annotation
varcount_sums_DN = "varcount_sums_DN"
#varcount_sums_СH = "varcount_sums_СH"

# This has to go, not used anywhere
gene_variant_mu_mask = "mus_DN"

# Mutational target output identifier
muttargs_DN_ID = "muttargs_DN"
muttargs_CH_ID = "muttargs_CH"

# Identifier for the QC file with per-proband variant count distributions
varcount_DN_mask = "variant_count_DN"
varcount_CH_mask = "variant_count_CH"
varcount_mask = "variant_count"

# Mutational targets per proband
muttargs_list_DN_ID = "muttargs_list_DN_ID"

mutnum_prod_CH = "mutnum_prod_CH"

mutnum_prod_dist_CH = "mutnum_prod_dist_CH"

DN_result = "DN_result"
CH_result = "CH_result"
CH_IND_result = "CH_IND_result"

# Precision for the p-value calculated with an infinite sum
pval_precision = 1e-6

# Binomial probability of the false diagnosis rate being higher 
#	than the given value
false_diag_rate = 0.05

# Maximal number of values calculated in the 'infinite' sum
maxIHval = 1000

# Threshold for the Irwin-Hall parameter above which the
#	Irwin-Hall distribution is approximated with Normal
IH_norm_approx_thr = 10


# Written on 01.23.2024 by Mikhail Moldovan, HMS DBMI