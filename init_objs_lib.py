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
# Gene_inst_dict: variant annotation -> ENSEMBL ID -> Gene instance (as specified in stat_lib)
def parse_variant_scores_file(variant_annot,
							  Gene_inst_dict,
							  score_thr):

	# total_mu is a normalizing value for mutational targets
	# After the parsing, all mutational targets are normalized to sum up to unit
	total_mu = 0
	# Input files are specified in cfg.py
	# variant_annot[0]: 'C'/'I' for coding or intronic variant types
	# variant_annot[1]: 'S'/'I', for SNP or indel variant types
	infile = cfg.variant_scores_files[variant_annot[0]][variant_annot[1]]
	with gzip.open(infile, 'rt') as inh:
		print(f"Parsing per-gene score-mutation rate distributions from: {infile}")
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
	Gene_inst_dict[variant_annot][current_ID] = sl.Gene(current_score_list,
														current_mu_list)

	total_mu += Gene_inst_dict[variant_annot][current_ID].mut_targ_arr[0]

	# Normalization of mutational targets
	for ENS_ID, Gene_obj in Gene_inst_dict[variant_annot].items():
		Gene_inst_dict[variant_annot][ENS_ID].normalize_mut_targs(total_mu)



# Function that parses precomputed mutational target files for CADD and SpliceAI score bins
# Produces: Gene_inst_dict: variant annotation -> ENSEMBL ID -> Gene instance
# total_mu_dict: variant annotation -> total mutational target of that annotation
def parse_variant_scores_files(score_thr_dict,
							   consequence_list,
			  				   suppress_indels):
	# Initialization of Gene_inst_dict
	# Gene_inst_dict: 
	#	variant annotation -> ENSEMBL ID -> Gene instance (as specified in stat_lib)
	Gene_inst_dict = {variant_annot : {} for variant_annot in variant_type_iter()}

	for variant_annot in variant_type_iter():
		# Not loading variant types that were not specified
		if variant_annot[0] not in consequence_list:
			continue

		if suppress_indels and variant_annot[1] == 'I':
			continue

		score_thr = score_thr_dict[variant_annot[0]]
		# A lot of code in this tool has the same structure:
		# Function that parses a list of files -> 
		#	function parsing a single file ->
		#	(sometimes) some functions processing a single line
		parse_variant_scores_file(variant_annot, 
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


# Makes a dictionary ensembl_id_dict: ENS_ID -> Gene ID
# Path to the input file is specified in cfg.py
def make_ENS2GeneID_dict(forward=True):
	# forward=True: ENS_ID -> Gene ID
	# forward=False: Gene ID -> ENS_ID
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
def load_s_het_bins():
	shet_dict = {}
	with gzip.open(cfg.shet_table, 'rt') as inh:
		for s in inh:
			if s.startswith('#'):
				continue
			s = s.strip().split()
			shet_dict[s[2]] = eval(s[11])
	# s_het_dict: ENSEMBL_ID -> s_het enrichment
	return shet_dict


# Written on 01.23.2024 by Mikhail Moldovan, HMS DBMI