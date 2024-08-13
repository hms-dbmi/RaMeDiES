"""
init_objs_lib

Additional library for RaMeDiES scripts containing mostly functions
for the parsing of default inputs.

Requirements:
cfg.py library
stat_lib.py library (specifically the Gene class)
python gzip library
"""

import cfg
import gzip
import stat_lib as sl

# Iterator over variant types
# variant type is a tuple of two values:
# 1. 'C'/'I' for coding or intronic variant
# 2. 'S'/'I', for SNP or indel variant
def variant_type_iter():
	# va: variant annotation (variant type)
	for va in cfg.var_annot_list:
		yield va

# Initiates variant count dictionary (varcount_dict)
# variant type (as specified in cfg.py and above) -> variant count
def init_varcount_dict():
	return {va : 0 for va in variant_type_iter()}


# Function that parses a precomputed mutational target file for CADD or SpliceAI score bins
# gene_instances: variant annotation -> ENSEMBL ID -> Gene instance (as specified in stat_lib)
def parse_variant_scores_file(variant_annot,
							  var_score_file_name,
							  Gene_inst_dict,
							  score_thr):

	# total_mu is a normalizing value for mutational targets
	# After the parsing, all mutational targets are normalized to sum up to unit
	total_mu = 0
	with gzip.open(var_score_file_name, 'rt') as inh:
		print(f"Parsing per-gene score-mutation rate distributions from: {var_score_file_name}")
		current_ID = None
		for inh_str in inh:
			if inh_str.startswith('#'):
				continue

			if inh_str.startswith("ENS"):
				if current_ID and current_score_list:
					# Initialization of Gene objects (specified in stat_lib)
					Gene_inst_dict[variant_annot][current_ID] = sl.Gene(current_score_list,
																		current_mu_list)

					# total_mu (mutational target of the genome) is calculated by summing 
					# 	mutational targets of lowest per-gene scores
					total_mu += Gene_inst_dict[variant_annot][current_ID].mut_targ_arr[0]
				current_ID = inh_str.strip()
				# current_mu_list: list of mutational targets
				# current_score_list: list of scores
				current_mu_list = []
				current_score_list = []
				continue

			if not current_ID:
				raise AssertionError(f"Wrong variant_scores_file format: {infile}")

			inh_str = inh_str.strip().split()
			score = eval(inh_str[0])
			# Not counting variants below a certain score threshold
			if score < score_thr:
				continue
			current_score_list.append(score)
			current_mu_list.append(eval(inh_str[1]))

	# Processing of the final gene
	if current_score_list:
		Gene_inst_dict[variant_annot][current_ID] = sl.Gene(current_score_list,
															current_mu_list)

		total_mu += Gene_inst_dict[variant_annot][current_ID].mut_targ_arr[0]

	# Normalization of mutational targets
	for ENS_ID, Gene_obj in Gene_inst_dict[variant_annot].items():
		Gene_inst_dict[variant_annot][ENS_ID].normalize_mut_targs(total_mu)



# Function that parses precomputed mutational target files for CADD and SpliceAI score bins
# Produces: gene_instances: variant annotation -> ENSEMBL ID -> Gene instance
# total_mu_dict: variant annotation -> total mutational target of that annotation
def parse_variant_scores_files(score_thr_dict,
							   consequence_list,
							   coding_score,
			  				   suppress_indels):
	# Initialization of gene_instances
	# gene_instances:
	#	variant annotation -> ENSEMBL ID -> Gene instance (as specified in stat_lib)
	Gene_inst_dict = {variant_annot : {} for variant_annot in variant_type_iter()}

	for variant_annot in variant_type_iter():
		# Not loading variant types that were not specified
		if variant_annot[0] not in consequence_list:
			continue

		if suppress_indels and variant_annot[1] == 'I':
			continue

		# Separate processing of Coding InDels
		if variant_annot == ('C', 'I'):
			score_thr = score_thr_dict["CInd"]
		else:
			score_thr = score_thr_dict[variant_annot[0]]

		# Input files are specified in cfg.py
		# variant_annot[0]: 'C'/'I' for coding or intronic variant types
		# variant_annot[1]: 'S'/'I', for SNP or indel variant types
		# coding_score: name of the score for coding SNPs provided under 
		#	--coding_score input parameter
		if variant_annot == ('C', 'S'):
			var_score_file_name = cfg.variant_scores_files['C']['S'][coding_score]
		else:
			var_score_file_name = cfg.variant_scores_files[variant_annot[0]][variant_annot[1]]
		# A lot of code in this tool has the same structure:
		# Function that parses a list of files -> 
		#	function parsing a single file ->
		#	(sometimes) some functions processing a single line
		parse_variant_scores_file(variant_annot,
								  var_score_file_name, 
								  Gene_inst_dict, 
								  score_thr)

	return Gene_inst_dict


# Name "pseudogene_dict" may be misleading
# This function parse a list of pseudogenes, RNA genes and overlapping genes
# The path to this list is specified in cfg.py
# Returns dictionary pseudogene_dict_ens: ENSEMBL ID -> True
# Python dictionary allows for a faster search
def make_pseudogene_dict():
	pseudogene_dict_ens = {}
	print(f"Parsing the list of discarded genes from: {cfg.pseudogenes}")
	with gzip.open(cfg.pseudogenes, 'rt') as inh:
		for s in inh:
			s = s.strip()
			pseudogene_dict_ens[s] = True
	return pseudogene_dict_ens


# Makes a dictionary ensembl_id_dict: ensembl_gene_id -> Gene ID
# Path to the input file is specified in cfg.py
def make_ENS2GeneID_dict(forward=True):
	# forward=True: ensembl_gene_id -> Gene ID
	# forward=False: Gene ID -> ensembl_gene_id
	ensembl_id_dict = {}
	with gzip.open(cfg.ens2gene, 'rt') as ens_handle:
		for id_line in ens_handle:
			id_line = id_line.strip().split()
			for ens_id in id_line[0].split(','):
				if forward:
					ensembl_id_dict[ens_id] = id_line[1]
				else:
					ensembl_id_dict[id_line[1]] = ens_id
	return ensembl_id_dict

# Loads s_het bin enrichment values from a file specified in cfg.py
# Enrichment values are the proportions of exclusively dominant OMIM genes
#	in s_het deciles
def load_gene_score_bins(scoreID = "s_het_R"):
	gene_score_dict = {}
	with gzip.open(cfg.gene_score_table, 'rt') as inh:
		for s in inh:
			s = s.strip().split()
			if s[0] == 'ensembl_gene_id':
				if scoreID not in s:
					raise AssertionError(f"ERROR: Gene score ID {scoreID} misspecified")
				score_index = s.index(scoreID)
				continue

			gene_score_dict[s[0]] = eval(s[score_index])
	# gene_score_dict: ENSEMBL_ID -> score bin enrichment for dominant disease-
	#	causing genes
	return gene_score_dict

# Write parameters used in the run to an intermediate file
def write_run_info(file_obj, args_obj):
	"""
	:param file_obj: python file object corresponding to an intermediate output file
	:args_obj: argparse arguments object with run specifications supplied by the user
	"""
	file_obj.write(f"### Information about RaMeDiES run ID {args_obj.o}\n")
	file_obj.write("### Do not discard this header when sharing intermediate outputs\n")
	file_obj.write(f"### VARIANT ANNOTATIONS: {args_obj.variant_annots}\n")
	file_obj.write(f"### INDEL SUPPRESSION: {args_obj.suppress_indels}\n")
	file_obj.write(f"### ONLY MISSENSE CODING SNPS: {args_obj.missense_run}\n")
	file_obj.write(f"### MAF CUTOFF: {args_obj.MAF}\n")
	if 'C' in args_obj.variant_annots:
		file_obj.write(f"### CODING SNP SCORE: {args_obj.coding_score}\n")
		file_obj.write(f"### CODING SNP SCORE THRESHOLD: {args_obj.C_thr}\n")
	if 'I' in args_obj.variant_annots:
		file_obj.write(f"### INTRONIC SCORE: SpliceAI\n")
		file_obj.write(f"### INTRONIC SCORE THRESHOLD: {args_obj.SAI_thr}\n")
	if not args_obj.suppress_indels:
		file_obj.write(f"### CODING INDEL SCORE: CADD\n")
		file_obj.write(f"### CODING INDEL SCORE THRESHOLD: {args_obj.CI_thr}\n")
	


# Written on 01.23.2024 by Mikhail Moldovan, HMS DBMI