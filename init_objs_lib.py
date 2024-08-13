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


def variant_type_iter():
	"""
	variant_annotation_type is a tuple with two values:
	1. C/I for coding/intronic variant
	2. S/I for SNV/indel variant
	:return: iterator over the variant annotation types present in the var_annot_list
	"""

	for variant_annotation_type in cfg.var_annot_list:
		yield variant_annotation_type


def init_varcount_dict():
	"""
	variant_annotation_type is a tuple with two values:
	1. C/I for coding/intronic variant
	2. S/I for SNV/indel variant
	:return: an empty variant count dictionary, {(C,S) -> 0, (C,I) -> 0, (I,S) -> 0, (I,I) -> 0}
	"""

	return {variant_annotation_type: 0 for variant_annotation_type in variant_type_iter()}


def parse_variant_scores_file(variant_annot, var_score_file_name, Gene_inst_dict, score_thr):
	"""
	Parse a specific precomputed mutational target file. Called by "parse_variant_scores_files"

	:param variant_annot: string corresponding to the variant annotation type (e.g., C for coding, S for SNVs)
	:param var_score_file_name: input file name (e.g., SpliceAI_IS.txt.gz)
	:param Gene_inst_dict: dictionary of variant_annotation_type -> Ensembl gene ID -> Gene instance (see stat_lin)
	:param score_thr: the minimal score threshold for variants considered
	:return: None, but update the Gene_inst_dict variable that was passed in
	"""

	total_mu = 0  # used to normalize the mutational targets so that they sum to 1

	with gzip.open(var_score_file_name, 'rt') as in_handle:
		print(f"Parsing per-gene score-mutation rate distributions from: {var_score_file_name}")

		current_ensg_id = None
		current_score_list = None  # list of deleteriousness scores
		current_mu_list = None  # list of corresponding mutational targets for those deleteriousness scores

		for in_handle_line in in_handle:
			if in_handle_line.startswith('#'):
				continue

			if in_handle_line.startswith("ENS"):

				if current_ensg_id and current_score_list:

					# store the values that were collected for the previous gene ID
					Gene_inst_dict[variant_annot][current_ensg_id] = sl.Gene(current_score_list, current_mu_list)

					# keep sum of total mu for this gene
					total_mu += Gene_inst_dict[variant_annot][current_ensg_id].mut_targ_arr[0]

				current_ensg_id = in_handle_line.strip()  # set the new Ensembl gene ID
				current_score_list = []  # clear out the score list
				current_mu_list = []  # clear out the corresponding mu list

				continue

			# we will never reach this point IF we have a properly formatted precomputed mutational target file
			if not current_ensg_id:
				raise AssertionError(f"Wrong variant_scores_file format: {var_score_file_name}")

			in_handle_line = in_handle_line.strip().split()  # score, mutational target

			score = eval(in_handle_line[0])
			if score < score_thr:  # skip variants below a certain score threshold
				continue

			current_score_list.append(score)
			current_mu_list.append(eval(in_handle_line[1]))

	# Process the final gene
	if current_score_list:
		Gene_inst_dict[variant_annot][current_ensg_id] = sl.Gene(current_score_list, current_mu_list)
		total_mu += Gene_inst_dict[variant_annot][current_ensg_id].mut_targ_arr[0]

	# Normalize the mutational targets so that they sum to 1 across all genes
	for ensembl_gene_id, Gene_obj in Gene_inst_dict[variant_annot].items():
		Gene_inst_dict[variant_annot][ensembl_gene_id].normalize_mut_targs(total_mu)


def parse_variant_scores_files(score_thr_dict, consequence_list, coding_score, suppress_indels):
	"""
	Parse the relevant set of precomputed mutational target files for score bins (e.g., CADD, REVEL, SpliceAI...)

	:param score_thr_dict: dictionary of variant_type (e.g., C=coding SNVs, CInd=coding indels) -> score threshold
	:param consequence_list: set of {C, I} indicating if coding and/or intronic regions are to be considered
	:param coding_score: name of the score for coding SNVs, provided by user using the --coding_score parameter
	:param suppress_indels: boolean indicating if indels are to be ignored (True)
	:return: dictionary of variant_annotation_type -> Ensembl gene ID -> Gene instance (see stat_lib)
	"""

	# initialize an empty dictionary of variant_annotation_type -> Ensembl gene ID -> Gene instance
	Gene_inst_dict = {variant_annot : {} for variant_annot in variant_type_iter()}

	for variant_annot in variant_type_iter():  # e.g., [('C', 'S'), ('C', 'I'), ('I', 'S'), ('I', 'I')]

		if variant_annot[0] not in consequence_list:  # skip loading any variants from regions that are excluded
			continue

		if suppress_indels and variant_annot[1] == 'I':  # skip loading indels if they are excluded
			continue

		if variant_annot == ('C', 'I'):  # process coding indels separately
			score_thr = score_thr_dict["CInd"]
		else:
			score_thr = score_thr_dict[variant_annot[0]]  # C=coding_SNVs, I=indels (same threshold for SNVs/indels)

		if variant_annot == ('C', 'S'):  # coding SNVs may have unique deleteriousness scores assigned
			var_score_file_name = cfg.variant_scores_files['C']['S'][coding_score]
		else:
			var_score_file_name = cfg.variant_scores_files[variant_annot[0]][variant_annot[1]]

		# parse the correct file for this variant annotation type and update the Gene_inst_dict accordingly
		parse_variant_scores_file(variant_annot, var_score_file_name,  Gene_inst_dict,  score_thr)

	return Gene_inst_dict


def make_pseudogene_dict():
	"""
	Create a list of Ensembl gene IDs that correspond to pseudogenes, RNA genes, and overlapping genes that
	will be EXCLUDED from our analysis!

	NOTE: the file path to the precomputed pseudogenes list is hardcoded in cfg

	:return: dictionary Ensembl gene ID -> boolean indicating skip status
	"""

	pseudogene_dict_ens = {}

	print(f"Parsing the list of discarded genes from: {cfg.pseudogenes}")
	with gzip.open(cfg.pseudogenes, 'rt') as in_handle:
		for ensembl_gene_id in in_handle:
			ensembl_gene_id = ensembl_gene_id.strip()
			pseudogene_dict_ens[ensembl_gene_id] = True

	return pseudogene_dict_ens


def make_ENS2GeneID_dict(forward=True):
	"""
	NOTE: path to precomputed mapping file is hardcoded in cfg

	:param forward: boolean indicating dictionary direction
	                True = ensembl_gene_id -> gene_name; False = gene_name -> ensembl_gene_id
	:return: dictionary mapping ensembl gene IDs to HGNC gene names
	"""

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


def load_gene_score_bins(constraint_score_id="s_het_R"):
	"""
	Briefly, all genes are sorted in ascending order according to their gene constraint value (specified by
	"constraint_score_id"). The sorted genes are split into deciles, and within each decile bin, the proportion of
	exclusively dominant OMIM genes is computed. This value is used to compute weighted Q-values which are then
	fed into an FDR procedure in RaMeDiES-DN.

	NOTE: path to precomputed gene constraint score bins is hardcoded in cfg

	:param constraint_score_id: type of gene constraint score to consider for weighted FDR procedure
	:return: mapping from Ensembl gene ID -> proportion of dominant OMIM genes in the score bin that gene falls into
	"""

	gene_score_dict = {}
	with gzip.open(cfg.gene_score_table, 'rt') as inh:
		for s in inh:
			s = s.strip().split()
			if s[0] == 'ensembl_gene_id':
				if constraint_score_id not in s:
					raise AssertionError(f"ERROR: Gene score ID {constraint_score_id} misspecified")
				score_index = s.index(constraint_score_id)
				continue

			gene_score_dict[s[0]] = eval(s[score_index])

	return gene_score_dict


def write_run_info(file_obj, args_obj):
	"""
	Write parameters used in the run to an intermediate file

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
		file_obj.write(f"### CODING SNP SCORE THRESHOLD: {args_obj.coding_snv_thr}\n")
	if 'I' in args_obj.variant_annots:
		file_obj.write(f"### INTRONIC SCORE: SpliceAI\n")
		file_obj.write(f"### INTRONIC SCORE THRESHOLD: {args_obj.SAI_thr}\n")
	if not args_obj.suppress_indels:
		file_obj.write(f"### CODING INDEL SCORE: CADD\n")
		file_obj.write(f"### CODING INDEL SCORE THRESHOLD: {args_obj.coding_indel_thr}\n")
	


# Written on 01.23.2024 by Mikhail Moldovan, HMS DBMI