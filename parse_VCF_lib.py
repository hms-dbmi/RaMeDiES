"""
parse_VCF_lib.py

Contains functions that parse a given VCF
Returns a set of updated Gene objects (specified in statlib.py) that are
	added to the [supplied] Gene_obj_dict dictionary containing information 
	about the cohort.

Dependencies:
1. cfg.py
2. stat_lib.py
3. python gzip library
4. python sys library
5. numpy

Main function is parse_VCF, which processes a single VCF file.
The rest are functions processing particular variants in the input VCF.
"""

import cfg
import stat_lib as sl
import init_objs_lib as iol

import gzip
import sys
import numpy as np

# Functions that checks that some needed column names are present in the header of
#	the processed VCF.
# The column names are specified in cfg.vcf_format_dict dictionary
# With different run specifications, some columns may be not needed.
# For instance, if processing of intronic regions is not specified,
#	the existence of SpliceAI score columns will be ignored.
def get_missing_terms(header, 
					  de_novo_bool, 
					  consequence_list, 
					  MAF_thr, 
					  no_qual_track_bool):
	
	# header: dictionary with column name codes (specified in cfg.vcf_format_dict)
	#	leading to column numbers in a tab-delimited VCF
	# de_novo_bool: boolean value specifying de novo run
	# consequence_list: list/string of single-character vriant position annotations
	#	'C' for coding, 'I' for intronic
	# MAF_thr: MAF threshold value, if not specified, the filter will not be imposed
	# no_qual_track_bool: boolean value for the case of no QC track in provided
	# VCFs

	# Initial missing terms list that is then filtered down according to
	#	run specifications
	missing_terms = []
	for VCF_field in cfg.vcf_format_dict.values():
		if not header.get(VCF_field):
			missing_terms.append(VCF_field)

	# Second missing terms list that is returned
	missing_terms_2 = []
	for VCF_field in missing_terms:
		# If coding variants are not processed, the CADD (or generally coding) 
		#	scores are not needed
		if 'C' not in consequence_list and cfg.rev_vcf_format_dict[VCF_field] == "CADD":
			continue
		# If intronic variants are not processed, the SpliceAI scores are not needed
		elif 'I' not in consequence_list and \
			cfg.rev_vcf_format_dict[VCF_field].startswith("SAI"):

			continue

		# For run on de novo VCFs, inheritance does not have to be specified
		if de_novo_bool and cfg.rev_vcf_format_dict[VCF_field] == "inherited_from":
			continue
		# If MAF threshold is not specified, the MAF column is not needed
		if not MAF_thr and cfg.rev_vcf_format_dict[VCF_field] == "MAF":
			continue
		# If no QC tracks are present in the VCFs, quality tracks will not be
		#	taken into account
		if no_qual_track_bool and cfg.rev_vcf_format_dict[VCF_field] == "qual_track":
			continue

		missing_terms_2.append(VCF_field)

	return ','.join(missing_terms_2)


# Function that infers coding/intronic annotation type from VEP consequences
def get_var_type1(var_annots, consequence_list):
	# var_type1: coding/intronic
	# Annotations may be separated by either comma (,) or an ampersand (&)
	var_annots = var_annots.replace('&', ',')
	conss = var_annots.split(',')
	# consequence codes are specified in cfg.VEP_cons_dict
	# Will crash the script to if VEP consequences are misspecified
	conss_codes = set([cfg.VEP_cons_dict[cons] for cons in conss])
	# Single annotation is assigned hierarchically:
	#	First, if coding consequences are present, variant is deemed coding
	#	Second, if no coding consequences are present and there are intronic 
	#	annotations, variant is deemed intronic
	if 'C' in conss_codes:
		var_type1 = 'C'
	elif 'I' in conss_codes:
		var_type1 = 'I'
	else:
		return None

	return var_type1

# Function that infers SNP/indel annotation
# Returns None if the length requirements for indels are not satisfied
def get_var_type2(ref_al, alt_al, var_type1, suppress_indels_bool):
	# var_type2: SNP/indel
	# ref_al/alt_al: reference and alternative alleles as specified in VCFs
	# var_type1 (coding/intronic) influences the length of accepted indels
	# suppress_indels_bool: under indel suppression, indels are not accepted

	# Mutation 'length': 0 for SNPs, negative for deletions, positive
	#	for insertions
	l_mut = len(ref_al) - len(alt_al)
	if l_mut == 0:
		return 'S'

	if l_mut != 0 and suppress_indels_bool:
		return None

	# Coding indels of the length of up to 10 are accepted
	if var_type1 == "C" and np.abs(l_mut) > 10:
		return None

	# Intronic deletions up to length 4 and 1nt insertions
	#	are accepted
	elif var_type1 == "I" and (l_mut > 4 or l_mut < -1):
		return None

	return 'I'

# Function that parses score values out of a VCF line
#	Returns None if the score is below the specified threshold 
def get_score_val(var_l, header, var_type1, score_thr_dict):
	# var_l: VCF line corresponding to a single variant
	# header: dictionary with column name codes (specified in cfg.vcf_format_dict)
	#	leading to column numbers in a tab-delimited VCF
	# var_type1: coding/intronic
	# score_thr_dict: dictionary with specified score thresholds:
	#	var_type1 (C/I for coding/intronic) -> score threshold
	if var_type1 == 'C':
		score = eval(var_l[header[cfg.vcf_format_dict["CADD"]]])
	else:
		# Maximal SpliceAI score over putative donor/acceptor gains/losses
		score = np.max([eval(var_l[header[cfg.vcf_format_dict["SAI_AG"]]]),
						eval(var_l[header[cfg.vcf_format_dict["SAI_AL"]]]),
						eval(var_l[header[cfg.vcf_format_dict["SAI_DG"]]]),
						eval(var_l[header[cfg.vcf_format_dict["SAI_DL"]]])])

	# Filter for scores below specified thresholds
	if score < score_thr_dict[var_type1]:
		return None
	return score


# Function that infers inheritence pattern of a given variant
#	from a line in a VCF
def get_inheritance(var_l, header, de_novo_bool, drop_HR = cfg.drop_HR):
	# var_l: VCF line corresponding to a single variant
	# header: dictionary with column name codes (specified in cfg.vcf_format_dict)
	#	leading to column numbers in a tab-delimited VCF
	# de_novo_bool: boolean value specifying de novo run
	# drop_HR: boolean value specifying homozyhous variant processing
	#	default False (not processed), specified in cfg.drop_HR

	# For a de novo run, the inherited_from column may be abscent in an 
	#	de novo-only VCF
	if de_novo_bool and not header.get(cfg.vcf_format_dict["inherited_from"]):
		return "DN"

	inh_info = var_l[header[cfg.vcf_format_dict["inherited_from"]]]

	# by the format specification (see manual), homozygous variant inheritance
	# consists of maternal/paternal inheritance keywords separated by a comma
	inh_info = inh_info.split(',')

	# Homozygous recessive variants are disregarded at this point
	if drop_HR and len(inh_info) > 1:
		return None

	# In-program inheritance identifier specified in cfg.inherited_from_dict
	inh_ID = cfg.inherited_from_dict.get(inh_info[0])
	if not inh_info:
		w_str = f"WARNING: Incorrect inheritance specification: {inh_info[0]}\n"
		sys.stderr.write(w_str)

	# Check if a variant is a de novo and if that agrees with the run specifications
	if de_novo_bool and inh_ID != "DN": 
		return None

	if not de_novo_bool and inh_ID == "DN": 
		return None

	return inh_ID


# Function that processes a quality track that may be present in a VCF.
def qual_track(var_l, header, no_qual_track_bool):
	# var_l: VCF line corresponding to a single variant
	# header: dictionary with column name codes (specified in cfg.vcf_format_dict)
	#	leading to column numbers in a tab-delimited VCF
	# no_qual_track_bool: boolean value for the abscence of quality track.
	#	is specified by the user as a parameter of one of master scripts

	if no_qual_track_bool:
		return True
	
	# quality track values are specified in cfg.vcf_format_dict.
	#	some value, followed by a comma, followed by high/low by default
	qual_val = var_l[header[cfg.vcf_format_dict["qual_track"]]].split(',')[-1]
	# Misspecification of quality values will produce a warning
	#	and the variant will be discarded
	qual_ID = cfg.qual_value_dict.get(qual_val)
	if qual_ID == None:
		sys.stderr.write(f"WARNING: Quality track misspecified: {qual_val}\n")
		return False
	return cfg.qual_value_dict[qual_val]


# Main function of the library. Processes a given VCF, returns:
#	1. Gene object updated with the de novo/CH variants from this VCF
#	2. varcount_dict: dictionary containing variant counts stratified by 
#		annotations (coding/intronic, SNP/indel) and inheritance patterns
# For the sake of computational speed, the filters are imposed sequentially
def parse_VCF(VCF_name,
			  Gene_inst_dict, 
			  gzip_bool,
			  MAF_thr,
			  score_thr_dict, 
			  pseudogene_dict,
			  consequence_list,
			  suppress_indels_bool,
			  de_novo_bool,
			  no_qual_track_bool):
	# VCF_name: path to the processed VCF
	# Gene_inst_dict: collection of Gene instances (class specified in stat_lib)
	#	stratified by annotations and ENSEMBL IDs
	# gzip_bool: boolean value indicating whether the VCF is gzip-compressed
	# MAF_thr: threshold for the minor allele frequency specified by user
	# score_thr_dict: collection of score thresholds specified by the user in 
	#	the master scripts
	# pseudogene_dict: list of discarded ENSEMBL IDs. Specified in 
	#	init_objs_lib/make_pseudogene_dict
	# consequence_list: variant consequence codes specified by the user.
	#	'C' for coding, 'I' for intronic
	# suppress_indels_bool: under indel suppression, indels are not accepted
	# de_novo_bool: boolean value specifying de novo run
	# no_qual_track_bool: boolean value for the abscence of quality track.


	# Selecting a proper file open function based on gzip_bool value
	(open_func, open_reg) = (gzip.open, 'rt') if gzip_bool else (open, 'r')

	# varcount_dict: inheritace code (DN/M/P) -> variant annotation -> variant count
	varcount_dict = {}

	# multiple_mut_dict: ENS_ID -> True. Used in the de novo processing to ensure
	#	only one de novo per gene per individual 
	multiple_mut_dict = {}
	with open_func(VCF_name, open_reg) as inh:
		# The first line in VCF starting with non-# is deemed a header
		h = None
		for var_l in inh:
			if var_l.startswith('#'):
				continue

			var_l = var_l.strip().split('\t')

			if not h:
				# h -> header: dictionary with column name codes 
				#	(specified in cfg.vcf_format_dict)
				h = {var_l[i]:i for i in range(len(var_l))}
				# Checking if the VCf format is consistent with the run setup
				missing_terms = get_missing_terms(h, 
												  de_novo_bool, 
												  consequence_list,
												  MAF_thr,
												  no_qual_track_bool)
				if missing_terms != '':
					sys.stderr.write(f"{VCF_name} Missing VCF fields: {missing_terms}\n")
					return None
				continue

#			print('\n', var_l)

			ENS_ID = var_l[h[cfg.vcf_format_dict["ENS_ID"]]]

			# Checking if the variant does not land within a discarded gene
			if pseudogene_dict.get(ENS_ID):
				continue

			# Infering standard parameters of a variant
			ref_al = var_l[h[cfg.vcf_format_dict["ref_al"]]]
			alt_al = var_l[h[cfg.vcf_format_dict["alt_al"]]]
			chrom = var_l[h[cfg.vcf_format_dict["chrom"]]]
#			print("1", chrom)

			# X-chromosome, Y-chromosome and mitochondrial variants are discarded
			if not chrom.isdigit():
				continue

			chrom = eval(chrom)
			position = eval(var_l[h[cfg.vcf_format_dict["position"]]])

			# Infering coding ('C')/intronic ('I') codes from VEP consequences
			var_annots = var_l[h[cfg.vcf_format_dict["var_annot"]]]
			var_type1 = get_var_type1(var_annots, consequence_list)
#			print("2", var_type1)

			# QC for consequence specification
			if not var_type1:
				continue

			# Infering whether the variant is an SNP ('S') or an indel ('I')
			var_type2 = get_var_type2(ref_al, alt_al, var_type1, suppress_indels_bool)
#			print("3", var_type2)

			# QC for indel length
			if not var_type2:
				continue

			# variant annotation as specified in cfg.var_annot_list
			var_annot = (var_type1, var_type2)

			# QC for information about this gene in the default inputs
			if not Gene_inst_dict[var_annot].get(ENS_ID):
				continue

			# Inference of CADD/SpliceAI scores
			score = get_score_val(var_l, h, var_type1, score_thr_dict)
#			print("4", score)

			# QC/thresholding of scores
			if score == None:
				continue

			# MAF filter
			if MAF_thr != -1:
				MAF = eval(var_l[header[cfg.vcf_format_dict["MAF"]]])
				if MAF > MAF_thr:
					continue

			# Inference of inheritance pattern
			inher = get_inheritance(var_l, h, de_novo_bool)
#			print("5", inher)

			# QC on inheritance
			if not inher:
				continue

			# Filter by the quality track
			if not qual_track(var_l, h, no_qual_track_bool):
				continue

			# Only the first of de novo variants occurring in a gene in a 
			#	patient will be processed
			if de_novo_bool and multiple_mut_dict.get(ENS_ID):
				sys.stderr.write(f"\nWARNING: multiple de novo variants in {VCF_name} gene {ENS_ID}\n")
				sys.stderr.write(f"...will include only the first observed variant\n")
				continue

			multiple_mut_dict[ENS_ID] = True

#			print("6", qual_track(var_l, h, no_qual_track_bool))
			
			# VCF name is used as individual ID throughout
			ind_id = VCF_name.split('/')[-1]

			# If all filters are passed, a Variant object (specified in stat_lib) is created
			Variant_obj = sl.Variant((chrom, position, ref_al, alt_al),
									  var_annot,
									  score,
									  inher,
									  ind_id)
			
			# ... and added to the vars list of a respective Gene object 
			#	(specified in stat_lib)
			Gene_inst_dict[var_annot][ENS_ID].vars.append(Variant_obj)

			# Updating the per-individual variant count dictionary
			# varcount_dict: inheritance code -> annotation code -> # variants
			if not varcount_dict.get(inher):
				varcount_dict[inher] = iol.init_varcount_dict()
			varcount_dict[inher][var_annot] += 1

	return varcount_dict


# Written on 01.23.2024 by Mikhail Moldovan, HMS DBMI