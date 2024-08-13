"""
parse_VCF_lib.py

Contains functions that parse a given tab-delimited variant input file
Returns a set of updated Gene objects (specified in stat_lib.py) that are added to the [supplied]
    Gene_obj_dict dictionary containing information about the cohort.

Dependencies:
1. cfg.py
2. stat_lib.py
3. python gzip library
4. python sys library
5. numpy

Main function is parse_variant_input.
"""

import cfg
import stat_lib
import init_objs_lib
import gzip
import sys
import numpy as np


# ------------------------------------------------------------------------------------------------------------
# Check that required columns are present in the header of the processed input variant files.
# Column names are specified in cfg.vcf_format_dict dictionary
# Some columns are optional depending on run specifications
# ------------------------------------------------------------------------------------------------------------

def get_missing_columns(header,
                        de_novo_bool,
                        suppress_indels_flag,
                        consequence_list,
                        maf_thr,
                        no_qual_track_bool):
    """
    :param header: dictionary with column name -> variable name
    :param de_novo_bool: True if running on de novos
    :param suppress_indels_flag: True if indels are to be ignored
    :param consequence_list: variant type annotations considered ('C' for coding, 'I' for intronic)
    :param maf_thr: no filter imposed if not specified
    :param no_qual_track_bool: False if quality control track is missing
    :return: comma-separated list of missing columns in input variant file
    """

    # total missing terms (filtered down according to run specifications)
    missing_terms = []
    for column_name in cfg.vcf_format_dict.values():
        if column_name not in header:
            missing_terms.append(column_name)

    # Reverse dictionary with input variant file format headers; note that some columns might be used more than once
    rev_vcf_format_dict = {}
    for k, v in cfg.vcf_format_dict.items():
        if v not in rev_vcf_format_dict:
            rev_vcf_format_dict[v] = set()
        rev_vcf_format_dict[v].add(k)

    missing_terms_2 = []

    # which columns do we care about?
    for column_name in missing_terms:

        # ignore coding scores if excluding coding variants
        if 'C' not in consequence_list and len([column_use for column_use in rev_vcf_format_dict[column_name]
                                                if column_use not in ["coding_snv_score", "coding_indel_score"]]) < 1:
            continue

        # ignore coding indel scores if indels are excluded:
        if suppress_indels_flag and len([column_use for column_use in rev_vcf_format_dict[column_name]
                                         if column_use != "coding_indel_score"]) < 1:
            continue

        # ignore SpliceAI scores if excluding intronic variants
        if 'I' not in consequence_list and len([column_use for column_use in rev_vcf_format_dict[column_name]
                                                if not column_use.startswith("SAI")]) < 1:
            continue

        # ignore inherited variant information when running in de novo mode
        if de_novo_bool and len([column_use for column_use in rev_vcf_format_dict[column_name]
                                 if column_use != "inherited_from"]) < 1:
            continue

        # ignore MAF column if no threshold is specified
        if maf_thr == -1 and len([column_use for column_use in rev_vcf_format_dict[column_name]
                                 if column_use != "MAF"]) < 1:
            continue

        # ignore quality Roulette track as specified
        if no_qual_track_bool and len([column_use for column_use in rev_vcf_format_dict[column_name]
                                 if column_use != "qual_track"]) < 1:
            continue

        missing_terms_2.append(column_name)

    return ','.join(missing_terms_2)


# ------------------------------------------------------------------------------------------------------------
# Infer if variant is coding or intronic from consequence column value
# ------------------------------------------------------------------------------------------------------------

def get_variant_loc_type(consequence_value,
                         missense_run_flag):
    """
    :param consequence_value: ,- or &-delimited values in the "consequence" column from the input variant file
    :param missense_run_flag: True means we are restricting to missense-only variants within exon regions
    :return: 'C' for coding variant, 'I' for intronic variant, None for neither
    """

    consequence_value = consequence_value.replace('&', ',').split(',')
    consequence_codes = set([cfg.VEP_cons_dict.get(cons, '') for cons in consequence_value])

    # Recall the following VEP consequence grouping:
    # 'U': UTR
    # 'I': intronic
    # 'S': synonymous
    # 'M': missense
    # 'C': coding
    # 'O': other

    if 'M' in consequence_codes:
        return 'C'
    elif 'C' in consequence_codes and not missense_run_flag:
        return 'C'
    elif 'I' in consequence_codes:
        return 'I'
    return None


# ------------------------------------------------------------------------------------------------------------
# Infer if variant is SNV or indel (within size limits) from ref and alt column values
# ------------------------------------------------------------------------------------------------------------

def get_variant_length_type(ref_al, alt_al, variant_loc_type, suppress_indels_bool):
    """
    :param ref_al: reference allele column value
    :param alt_al: alternate allele column value
    :param variant_loc_type: 'C' for coding and 'I' for intronic
    :param suppress_indels_bool: True if indels are to be ignored
    :return: 'S' for SNV, 'I' for indel, None for neither
    """

    variant_length = len(ref_al) - len(alt_al)

    if variant_length == 0:
        return 'S'
    if variant_length != 0 and suppress_indels_bool:
        return None
    if variant_loc_type == "C" and np.abs(variant_length) > 10:  # coding indels must be <=10bp
        return None
    elif variant_loc_type == "I" and (
            variant_length > 4 or variant_length < -1):  # intronic insertions must be <=4bp and deletions <=1 bp
        return None

    return 'I'


# ------------------------------------------------------------------------------------------------------------
# Get corresponding variant functionality score given the variant type
# ------------------------------------------------------------------------------------------------------------

def get_score_val(variant_line, header, variant_type, score_thr_dict):
    """
    :param variant_line: list of tab-delimited values from processed variant input file
    :param header: dictionary of column_name -> 0-index in file
    :param variant_type: pair of values: (0: 'C' for coding, 'I' for intronic, 1: 'S' for SNV', 'I' for indel)
    :param score_thr_dict: dictionary of variant type -> variant score threshold
    :return: float of the appropriate score
    """

    if variant_type[0] == 'C':
        if variant_type[1] == 'I':
            score = variant_line[header[cfg.vcf_format_dict["coding_indel_score"]]]
        else:
            score = variant_line[header[cfg.vcf_format_dict["coding_score"]]]
        if not score[-1].isdigit():
            return None
        else:
            score = eval(score)
    else:
        # Maximal SpliceAI score over putative donor/acceptor gains/losses
        score = np.max([eval(variant_line[header[cfg.vcf_format_dict["SAI_AG"]]]),
                        eval(variant_line[header[cfg.vcf_format_dict["SAI_AL"]]]),
                        eval(variant_line[header[cfg.vcf_format_dict["SAI_DG"]]]),
                        eval(variant_line[header[cfg.vcf_format_dict["SAI_DL"]]])])

    # Filter for scores below specified thresholds
    variant_type_str = cfg.var_type_to_var_str[variant_type]
    if score < score_thr_dict[variant_type_str]:
        return None
    return score


# ------------------------------------------------------------------------------------------------------------
# Get inheritance ID (M=maternally inherited, P=paternally inherited, DN=denovo)
# ------------------------------------------------------------------------------------------------------------

def get_inheritance(variant_line, header, de_novo_bool, drop_homozygous_recessive=cfg.drop_HR):
    """
    :param variant_line: list of tab-delimited values from processed variant input file
    :param header: dictionary of column_name -> 0-index in file
    :param de_novo_bool: True if de novos are being evaluated
    :param drop_homozygous_recessive: True if homozygous recessive variants should be processed
    :return: 'DN' if this is a de novo and None otherwise
    """

    if de_novo_bool and not header.get(cfg.vcf_format_dict["inherited_from"]):
        return "DN"

    inheritance_value = variant_line[header[cfg.vcf_format_dict["inherited_from"]]].split(',')

    # Homozygous recessive variants are disregarded at this point
    if drop_homozygous_recessive and len(inheritance_value) > 1:
        return None

    # In-program inheritance identifier specified in cfg.inherited_from_dict
    inheritance_id = cfg.inherited_from_dict.get(inheritance_value[0])
    if not inheritance_value:
        w_str = f"WARNING: Incorrect inheritance specification: {inheritance_value[0]}\n"
        sys.stderr.write(w_str)

    # only considering de novos (and this is inherited)?
    if de_novo_bool and inheritance_id != "DN":
        return None

    # only considering inherited (and this is a denovo) ?
    if not de_novo_bool and inheritance_id == "DN":
        return None

    return inheritance_id


# ------------------------------------------------------------------------------------------------------------
# Parse the Roulette quality track
# ------------------------------------------------------------------------------------------------------------

def qual_track(variant_line, header, no_qual_track_bool):
    """
    :param variant_line: list of tab-delimited values from processed variant input file
    :param header: dictionary of column_name -> 0-index in file
    :param no_qual_track_bool: True if Roulette quality track is missing
    :return: True for passing quality
    """

    if no_qual_track_bool:
        return True

    # e.g., 8.456e-10,high
    quality_value = variant_line[header[cfg.vcf_format_dict["qual_track"]]].split(',')[-1]
    quality_id = cfg.qual_value_dict.get(quality_value)

    if not quality_id:
        # sys.stderr.write('! Low quality variant dropped: ' + ' '.join(variant_line) + '\n')
        return False

    return quality_id


def is_gzipped(file_path):
    """
    Check if a file is gzipped by checking its file signature.
    """
    with open(file_path, 'rb') as f:
        # Read the first two bytes to check the file signature
        signature = f.read(2)
    return signature == b'\x1f\x8b'  # gzip file signature


# ------------------------------------------------------------------------------------------------------------
# MAIN FUNCTION
# for a given input variant file, return gene object with all variants and a dictionary of variant type -> counts
# ------------------------------------------------------------------------------------------------------------

def parse_variant_input(input_variant_file,
                        gene_instances,
                        maf_threshold,
                        score_threshold_dict,
                        pseudogene_dict,
                        consequence_list,
                        suppress_indels_flag,
                        de_novo_flag,
                        no_qual_track_flag,
                        missense_run_flag):
    """
    :param input_variant_file: full path to a processed tab-delimited variant file
    :param gene_instances: collection of Gene instances (class specified in stat_lib) stratified by annotations and
        ensembl IDs
    :param maf_threshold: user-specified minor allele frequency threshold
    :param score_threshold_dict: variant type -> score threshold imposed
    :param pseudogene_dict: ensembl IDs to be excluded, specified in init_objs_lib/make_pseudogene_dict
    :param consequence_list: variant consequences considered; C=coding, I=intronic
    :param suppress_indels_flag: True if indels are to be ignored
    :param de_novo_flag: True if we are processing de novos
    :param no_qual_track_flag: True if Roulette quality track is absent
    :param missense_run_flag: True if only missenses should be considered
    :return: dict of variant counts by annotation
    """

    # Selecting a proper file open function based on gzip_flag value
    (open_func, open_reg) = (gzip.open, 'rt') if is_gzipped(input_variant_file) else (open, 'r')

    varcount_dict = {}  # inheritace code (DN/M/P) -> variant annotation -> variant count

    multiple_mutation_dict = {}  # ensembl ID -> True if 2+ de novos in the same gene in the same individual

    with open_func(input_variant_file, open_reg) as inh:
        header = None
        for variant_line in inh:
            if variant_line.startswith('#'):
                continue

            variant_line = variant_line[:-1].split('\t')

            if not header:
                header = {column_name: column_index for column_index, column_name in enumerate(variant_line)}
                missing_terms = get_missing_columns(header,
                                                    de_novo_flag,
                                                    suppress_indels_flag,
                                                    consequence_list,
                                                    maf_threshold,
                                                    no_qual_track_flag)
                if missing_terms != '':
                    sys.stderr.write(f"{input_variant_file} Missing VCF fields: {missing_terms}\n")
                    return None
                continue

            ensembl_id = variant_line[header[cfg.vcf_format_dict["ensembl_gene_id"]]]

            # skip blacklisted pseudogenes
            if pseudogene_dict and ensembl_id in pseudogene_dict:
                continue

            # restrict to autosomes 1-22
            chrom = variant_line[header[cfg.vcf_format_dict["chrom"]]].replace('chr', '')
            if not chrom.isdigit():  # drop chromosomes X, Y, MT
                continue
            chrom = eval(chrom)

            # is variant a coding or intronic variant?
            variant_types = variant_line[header[cfg.vcf_format_dict["var_annot"]]]  # e.g., missense&nonsense

            # e.g., C for (missense) coding, I for intronic, None if neither
            var_type1 = get_variant_loc_type(variant_types, missense_run_flag)
            if not var_type1:
                continue

            # is variant an SNV or appropriately-sized indel?
            ref_al = variant_line[header[cfg.vcf_format_dict["ref_al"]]]
            alt_al = variant_line[header[cfg.vcf_format_dict["alt_al"]]]
            var_type2 = get_variant_length_type(ref_al, alt_al, var_type1, suppress_indels_flag)
            if not var_type2:
                continue

            # is the variant a type we are interested in?
            variant_type = (var_type1, var_type2)
            if gene_instances and ensembl_id not in gene_instances[variant_type]:
                continue

            # is the variant functionality score high enough?
            score = get_score_val(variant_line, header, variant_type, score_threshold_dict)
            if not score:
                continue

            # is the variant rare enough in the population?
            if maf_threshold != -1:
                maf = eval(variant_line[header[cfg.vcf_format_dict["MAF"]]])
                if maf > maf_threshold:
                    continue

            # is the inheritance pattern one we are interested in?
            inheritance = get_inheritance(variant_line, header, de_novo_flag)
            if not inheritance:
                continue

            # is this a high quality variant (if we care)?
            if not qual_track(variant_line, header, no_qual_track_flag):
                continue

            # select the FIRST de novo per gene per patient (if there are multiple)
            if de_novo_flag and multiple_mutation_dict.get(ensembl_id):
                sys.stderr.write(f"\nWARNING: multiple de novo variants in {input_variant_file} gene {ensembl_id}\n")
                sys.stderr.write(f"...will include only the first observed variant\n")
                continue
            multiple_mutation_dict[ensembl_id] = True

            # if we made it to this point, create a new variant object
            position = eval(variant_line[header[cfg.vcf_format_dict["position"]]])
            individual_id = input_variant_file.split('/')[-1]
            variant_object = stat_lib.Variant((chrom, position, ref_al, alt_al),
                                              variant_type,
                                              score,
                                              inheritance,
                                              individual_id)

            # and add it to the corresponding gene object
            gene_instances[variant_type][ensembl_id].vars.append(variant_object)

            # Updating the per-individual variant count dictionary
            # varcount_dict: inheritance code -> annotation code -> # variants
            if not varcount_dict.get(inheritance):
                varcount_dict[inheritance] = init_objs_lib.init_varcount_dict()
            varcount_dict[inheritance][variant_type] += 1

    return varcount_dict


# Written on 01.23.2024 by Mikhail Moldovan, HMS DBMI
