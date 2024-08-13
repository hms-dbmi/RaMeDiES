"""
stat_lib

Library containing core classes and functions used in 
ramediesDN, ramediesCH and ramedies_CH_IND.

##############
Dependencies:
1. numpy
2. os
3. scipy
4. init_objs_lib
5. cfg

##############
Classes:
Variant class: class with instances containing information about particular variants

Gene class: class containing score-mutational target distributions as well as 
    the variant information for each gene stratified by annotation.

CH_variant class: class storing information about CH variants in the form of pairs
    of Variant class instances.

VariantCollection class: instances store variant (Variant class) information
    stratified by individuals. Used in ramediesCH and ramedies_CH_IND algorithms

##############
Functions are subdivided into blocks. See the rest of the document for details.

"""

import numpy as np
import init_objs_lib
import cfg
from os import path, remove
from scipy.stats import norm, binom, poisson


# ==========================================================================================
# Universal Class Definitions
# ==========================================================================================

class Variant:
    """
    Variant contains relevant information about a particular variant.

    Attributes:
        mut_id_tuple (str, str, str, str): chromosome, 1-indexed position, reference allele, alternate allele.
        var_annot (char, char): variant type 'C' for coding or 'I' for intronic; 'S' for SNV or 'I' for indel.
        score (float): deleteriousness score of variant.
        inher (str): variant inheritance; M = from mom, D = from dad, DN = neither/de novo.
        proband_id (str): name of VCF file serving as a proband identifier.
    """

    def __init__(self, mut_id_tuple, var_annot, score, inher, proband_id=None):
        """
        Initializes Variant with provided variant attributes.

        Args:
            mut_id_tuple (str, str, str, str): chromosome, 1-indexed position, reference allele, alternate allele.
            var_annot (char, char): variant type 'C' for coding or 'I' for intronic; 'S' for SNV or 'I' for indel.
            score (float): deleteriousness score of variant.
            inher (str): variant inheritance; M = from mom, D = from dad, DN = neither/de novo.
            proband_id (str): name of VCF file serving as a proband identifier.
        """
        self.mut_id_tuple = mut_id_tuple
        self.var_annot = var_annot
        self.score = score
        self.inher = inher
        self.proband_id = proband_id

    def make_mut_id(self):
        """
        :return: (str) corresponding to the variant ID, formatted as chromosome:ref-allele_position_alt-allele
        """
        chrom = self.mut_id_tuple[0]
        ref = self.mut_id_tuple[2]
        pos = self.mut_id_tuple[1]
        alt = self.mut_id_tuple[3]
        mut_id = f"chr{chrom}:{ref}_{pos}_{alt}"
        return mut_id

    def print_info(self):
        """
        :return: (str) containing all information contained in the Variant class
        """
        mut_id = self.make_mut_id()
        va = f"{self.var_annot[0]}{self.var_annot[1]}"  # variant annotation e.g., CI for coding indel
        return f"{mut_id}|{va}|{self.score}|{self.inher}|{self.proband_id}"


class Gene:
    """
    Gene contains per-gene mutational target distributions by annotation class (coding/intronic and SNV/indel)
    Gene also stores observed variants as a list of Variant instances

    Attributes:
        score_arr [float]: array of sorted (ascending) values of unique deleteriousness scores present in gene
        mut_targ_arr [float]: array of corresponding mutational targets for each unique deleteriousness score
            computed as the total mutational target of all scores equal to or higher than that score
        gene_mu (float): the maximum mutational target of this gene (i.e., all variants included)
        gene_length (int): total number of unique deleteriousness scores in this gene, unrelated to actual gene length
        vars [Variant]: list of variants observed in this gene in a proband
        mu_list [float]: list of mutational targets of all variants of a specific type
    """
    def __init__(self, score_arr, mut_targ_arr):
        """
        Args:
            score_arr [float]: read from precomputed data/score_lists*gz files
            mut_targ_arr [float]: read from precomputed data/score_lists*gz files
        """
        self.score_arr, self.mut_targ_arr = zip(*sorted(zip(score_arr, mut_targ_arr)))
        self.mut_targ_arr = np.array(self.mut_targ_arr)
        self.gene_mu = self.mut_targ_arr[0]
        self.gene_length = len(self.mut_targ_arr)
        self.vars = []  # updated by the parse_VCF_lib/parse_variant_input function
        self.mu_list = []

    def normalize_mut_targs(self, mut_targ_norm=1.0):
        """
        Normalize mutational targets so that the mutational target of the entire genome will be 1.
        Normalizations are generated by parse_VCF_lib/parse_variant_scores_files function as total_mu_dict dictionary

        Args:
            mut_targ_norm (float): normalize mutational targets to sum to this value (default: 1.0)
        """
        self.mut_targ_arr = self.mut_targ_arr / mut_targ_norm
        self.gene_mu = self.mut_targ_arr[0]

    def get_mu_from_score(self, score):
        """
        :return: the mutational target in this gene corresponding to a given score (search binary search)

        Args:
            score (int): deleteriousness score to return a corresponding mutational target for
        """
        if self.score_arr[0] >= score:
            return self.mut_targ_arr[0]
        elif self.score_arr[-1] < score:
            return 0.0

        ind_low = 0
        ind_high = self.gene_length - 1
        ind_mid = (ind_high - ind_low) // 2

        while ind_high - ind_low > 1:
            if self.score_arr[ind_mid] == score:
                return self.mut_targ_arr[ind_mid]

            elif self.score_arr[ind_low] == score:
                return self.mut_targ_arr[ind_low]

            elif self.score_arr[ind_high] == score:
                return self.mut_targ_arr[ind_high]

            elif self.score_arr[ind_mid] < score:
                ind_low = ind_mid
                ind_mid = (ind_high + ind_low) // 2

            elif self.score_arr[ind_mid] > score:
                ind_high = ind_mid
                ind_mid = (ind_high + ind_low) // 2

        return self.mut_targ_arr[ind_low]

    def calculate_muttargs_denovo(self, sum_values=False):
        """
        :return: the mutational target for the y statistic for a single annotation (used only in the de novo regime)

        Args:
            sum_values (bool): return the sum of mutational targets rather than the list (default: False)
        """
        self.mu_list = []  # list of mutational targets which is printed into a log file
        for var in self.vars:
            mu = self.get_mu_from_score(var.score)
            self.mu_list.append((self.gene_mu - mu) / self.gene_mu)

        if sum_values:
            return np.sum(self.mu_list)

        return self.mu_list


# ==========================================================================================
# Compound Heterozygous Class Definitions (used by RaMeDiES-CH and RaMeDiES-IND)
# ==========================================================================================

class CH_variant:
    """
    CH_variant contains relevant information for a compound heterozygous variant pair

    Attributes:
        var_P (Variant): paternally inherited Variant
        var_M (Variant): maternally inherited variant
        ensembl_gene_id (str): Ensembl gene ID where the variant is found
        mu (float): maximum mutational target of maternally/paternally inherited variants
        mu2 (float): normalized squared mutational target corresponding to the compound heterozygous variant pair
    """
    def __init__(self, Variant_obj_tuple, ensembl_gene_id):
        """
        Args:
            Variant_obj_tuple ((Variant, Variant)): pair of paternally and maternally inherited variants in a gene
            ensembl_gene_id: Ensembl gene ID where the variants were found in trans
        """
        self.var_P = Variant_obj_tuple[0]
        self.var_M = Variant_obj_tuple[1]
        self.ENS_ID = ensembl_gene_id
        self.mu = None
        self.mu2 = None

    def muttarg(self, Gene_obj_dict, norm=None):
        """
        :return: mutational target of the compound heterozygous variant

        Args:
            Gene_obj_dict (dict): variant annotation (e.g., CI) -> Ensembl gene ID -> Gene (object)
            norm: normalization factor to scale squared mutational targets
        """
        # Gene_obj_dict: variant annotation -> ENSEMBL ID -> Gene object
        Gene_obj_P = Gene_obj_dict[self.var_P.var_annot][self.ENS_ID]
        Gene_obj_M = Gene_obj_dict[self.var_M.var_annot][self.ENS_ID]
        mu_P = Gene_obj_P.get_mu_from_score(self.var_P.score)
        mu_M = Gene_obj_M.get_mu_from_score(self.var_M.score)

        if not norm:
            norm = np.max([Gene_obj_P.gene_mu ** 2, Gene_obj_M.gene_mu ** 2])

        self.mu = np.max([mu_P, mu_M])  # unadjusted mutational target of compound heterozygous variant
        self.mu2 = self.mu ** 2 / norm  # compound heterozygous mutational target with respect to the gene

        return self.mu2

    def print_info(self):
        """
        :return: (str) corresponding to paternal variant & maternal variant information
        """
        info_P = self.var_P.print_info()
        info_M = self.var_M.print_info()
        return f"{info_P}&{info_M}"


class VariantCollection:
    """
    VariantCollection stores information about singular and compound heterozygous variants

    Attributes:
        P_vars ({Variant}): set of paternally inherited variants
        M_vars ({Variant}): set of maternally inherited Variants
        CH_list ([CH_variant]): list of all unique compound heterozygous variant pair objects in this gene
        ensembl_gene_id (str): Ensembl gene ID
        y (str): y statistic corresponding to the most surprising compound heterozygous variant mutational target
                 in this gene (scaled with respect to the genome)
    """
    def __init__(self, ensembl_gene_id, Gene_obj_dict):
        """
        Sort the variants present in Gene_obj_dict by inheritance

        Args:
            ensembl_gene_id (str): Ensembl gene ID
            Gene_obj_dict (dict): variant annotation (e.g., CI) -> Ensembl gene ID -> Gene (object)
        """
        self.P_vars = {}
        self.M_vars = {}
        for va in init_objs_lib.variant_type_iter():
            if not Gene_obj_dict[va].get(ensembl_gene_id):
                continue

            Gene_obj = Gene_obj_dict[va][ensembl_gene_id]
            for Variant_obj in Gene_obj.vars:

                if Variant_obj.inher == 'P':
                    if not self.P_vars.get(Variant_obj.proband_id):
                        self.P_vars[Variant_obj.proband_id] = []
                    self.P_vars[Variant_obj.proband_id].append(Variant_obj)

                elif Variant_obj.inher == 'M':
                    if not self.M_vars.get(Variant_obj.proband_id):
                        self.M_vars[Variant_obj.proband_id] = []
                    self.M_vars[Variant_obj.proband_id].append(Variant_obj)

        self.CH_list = []
        self.ENS_ID = ensembl_gene_id
        self.y = None

    def make_comphet_var_list(self, Gene_obj_dict):
        """
        :return: list of compound heterozygous variant with the smallest mutational target per gene in this proband

        Args:
            Gene_obj_dict (dict): variant annotation (e.g., CI) -> Ensembl gene ID -> Gene (object)
        """
        if self.P_vars == {} or self.M_vars == {}:
            return None

        # for each proband, find all compound heterozygous variant pairs (trans within a gene)
        per_proband_CHs = {}
        for proband_id, Variant_obj_list_P in self.P_vars.items():
            Variant_obj_list_M = self.M_vars.get(proband_id)
            if not Variant_obj_list_M:
                continue

            per_proband_CHs[proband_id] = {}

            for P_var in Variant_obj_list_P:
                for M_var in Variant_obj_list_M:
                    CH = CH_variant((P_var, M_var), self.ENS_ID)
                    mu2_CH = CH.muttarg(Gene_obj_dict)
                    per_proband_CHs[proband_id][mu2_CH] = CH

        for proband_id, CH_dict in per_proband_CHs.items():
            mu2_top = sorted(list(CH_dict.keys()))[0]
            self.CH_list.append(CH_dict[mu2_top])

    def calc_comphet_y(self):
        """
        :return: the cohort-wide y statistic
        """
        y = 0
        for CH_obj in self.CH_list:
            y += 1 - CH_obj.mu2

        self.y = y
        return y


# ==========================================================================================
# Write and load summary statistics files containing denovo VARIANT COUNT information
# ==========================================================================================

def write_by_annot_varcount(by_annot_varcount_dict, outfile_mask, args_obj):
    """
    :param by_annot_varcount_dict: inheritance type (e.g., M for maternal) -> variant type (e.g., CI for coding indels)
                                   -> variant count (number of variants of this type)
    :param outfile_mask: outfile file prefix, specified by user using the --o parameter
    :param args_obj: dictionary of parsed command-line arguments
    :return: None, write variant counts by variant type to outfile
    """

    outfile_name = f"{outfile_mask}_{cfg.varcount_sums_DN}.txt"

    with open(outfile_name, 'w') as outh:
        outh.write('# De novo mutation count by variant type across cohort\n')
        init_objs_lib.write_run_info(outh, args_obj)
        outh.write("inheritance\tvariant_type\tdenovo_mutation_count\n")
        for inher, var_count_dict in by_annot_varcount_dict.items():
            for var_annot, var_count in by_annot_varcount_dict[inher].items():
                va_str = f"{var_annot[0]}{var_annot[1]}"
                outh.write(f"{inher}\t{va_str}\t{var_count}\n")

    print(f"By-annotation variant counts written to: {outfile_name}")


def make_by_annot_varcount_dict(varcount_dict, outfile_mask, args_obj):
    """
    :param varcount_dict: individual ID (filename) -> inheritance type (e.g., M for maternal) ->
                          variant type (e.g., CI for coding indels) -> # variants of this type
    :param outfile_mask: outfile file prefix, specified by user using the --o parameter
    :param args_obj: dictionary of parsed command-line arguments
    :return: the "by_annot_varcount_dict" used as an argument to above function "write_by_annot_varcount"
    """

    by_annot_varcount_dict = {}
    for ind_id, ind_varcount in varcount_dict.items():
        for inher, inher_dict in ind_varcount.items():
            if not by_annot_varcount_dict.get(inher):
                by_annot_varcount_dict[inher] = init_objs_lib.init_varcount_dict()
            for va in cfg.var_annot_list:
                by_annot_varcount_dict[inher][va] += ind_varcount[inher][va]

    write_by_annot_varcount(by_annot_varcount_dict, outfile_mask, args_obj)
    return by_annot_varcount_dict


def load_varcount_from_file(by_annot_varcount_dict,
                            input_ID,
                            variant_annots,
                            suppress_indels_bool):
    """
    :param by_annot_varcount_dict: inheritance type (e.g., M for maternal) -> variant type (e.g., CI for coding indels)
                                   -> variant count (number of variants of this type)
    :param input_ID: summary statistics input file ID (initially produced by the "write_by_annot_varcount" function),
                     specified by user using the --M parameter
    :param variant_annots: variant annotation types included, specified by user using the --variant_annots parameter
    :param suppress_indels_bool: boolean indicating whether to include indels (False), specified by user
    :return: updated "by_annot_varcount_dict" used as an argument to above function "write_by_annot_varcount"
    """

    infile_name = f"{input_ID}_{cfg.varcount_sums_DN}.txt"

    if not path.isfile(infile_name):
        raise AssertionError(f"Variant count file {infile_name} does not exist")

    with open(infile_name, 'r') as inh:
        header = None
        for vc_str in inh:
            if vc_str.startswith('#'):
                continue
            if not header:
                header = vc_str.strip().split()
                continue

            vc_str = vc_str.strip().split()

            inher = vc_str[0]  # 1st column: inheritance code (cfg.inherited_from_dict)

            var_annot = (vc_str[1][0], vc_str[1][1])  # 2nd column: variant annotation (cfg.var_annot_list)
            if not var_annot[0] in variant_annots:
                continue
            if suppress_indels_bool and var_annot[1] == 'I':
                continue

            varcount = eval(vc_str[2])  # 3rd column: variant count

            if not by_annot_varcount_dict.get(inher):  # by_annot_varcount_dict update
                by_annot_varcount_dict[inher] = init_objs_lib.init_varcount_dict()
            by_annot_varcount_dict[inher][var_annot] += varcount

    print(f"Variant counts from {infile_name} loaded")
    return by_annot_varcount_dict


def load_by_annot_varcount_dict(input_IDs, variant_annots, suppress_indels_bool):
    """
    :param input_IDs: comma-separated list of summary statistics input file IDs (initially produced by the
                      "write_by_annot_varcount" function), specified by user using the --M parameter
    :param variant_annots: variant annotation types included, specified by user using the --variant_annots parameter
    :param suppress_indels_bool: boolean indicating whether to include indels (False), specified by user
    :return: updated "by_annot_varcount_dict" used as an argument to above function "write_by_annot_varcount"
    """

    input_IDs = input_IDs.split(',')
    by_annot_varcount_dict = {}
    for input_ID in input_IDs:
        by_annot_varcount_dict = load_varcount_from_file(by_annot_varcount_dict,
                                                         input_ID,
                                                         variant_annots,
                                                         suppress_indels_bool)

    return by_annot_varcount_dict


def write_varcount_dist(varcount_dict,
                        outfile_mask, 
                        input_directory='<unspecified>',
                        args_obj=None):
    """
    :param varcount_dict: individual ID (filename) -> inheritance type (e.g., M for maternal) ->
                          variant type (e.g., CI for coding indels) -> # variants of this type
    :param outfile_mask: outfile file prefix, specified by user using the --o parameter
    :param input_directory: directory containing all tab-delimited input variant files per individual
    :param args_obj: dictionary of parsed command-line arguments
    :return: None, write variant count distribution to file. NOTE: this is not used in any analysis, but may
             be useful for quality control analyses of variant calling across the cohort
    """

    varcount_dist = {}
    suffix = 'comphet'
    for file_name, inher_dict in varcount_dict.items():
        for inher, va_dict in inher_dict.items():
            if 'DN' in inher:
                suffix = 'denovo'
            for va, var_num in va_dict.items():
                if not varcount_dist.get(inher):
                    varcount_dist[inher] = {va: {} for va in cfg.var_annot_list}
                if not varcount_dist[inher][va].get(var_num):
                    varcount_dist[inher][va][var_num] = 0
                varcount_dist[inher][va][var_num] += 1

    outfile_name = f"{outfile_mask}_{suffix}_{cfg.varcount_mask}.txt"

    # Printing the variant count distribution to the output
    with open(f"{outfile_name}", 'w') as outh:
        outh.write('# Variant count distribution computed from input variant files in '+input_directory+'\n')
        init_objs_lib.write_run_info(outh, args_obj)
        outh.write("inheritance\tvariant_type\tvariant_count\tnumber_samples\n")
        for inher, va_dict in varcount_dist.items():
            for va, dist_dict in va_dict.items():
                va_str = f"{va[0]}{va[1]}"
                for var_num in sorted(dist_dict.keys()):
                    freq = dist_dict[var_num]
                    outh.write(f"{inher}\t{va_str}\t{var_num}\t{freq}\n")

    print(f"Variant count distributions written to: {outfile_name}")


# ==========================================================================================
# Write and load summary statistics files containing denovo MUTATIONAL TARGET information
# ==========================================================================================

def load_muttargs_from_file(Gene_inst_dict,
                            input_file_id,
                            variant_annots,
                            suppress_indels_bool):
    """
    Iteratively called by load_muttargs_from_filelist function below

    :param Gene_inst_dict: variant annotation type (e.g., CI) -> Ensembl Gene ID -> Gene (object)
    :param input_file_id: summary statistics input file ID (initially produced by the "write_muttargs" function),
                     specified by user using the --M parameter
    :param variant_annots: variant annotation types included, specified by user using the --variant_annots parameter
    :param suppress_indels_bool: boolean indicating whether to include indels (False), specified by user
    :return: updated (previously initialized) "Gene_inst_dict" used as an argument to below function "write_muttargs"
    """

    infile_name = f"{input_file_id}_{cfg.muttargs_list_DN_ID}.txt"
    if not path.isfile(infile_name):
        raise AssertionError(f"mutational target file {infile_name} does not exist")

    with open(infile_name, 'r') as inh:
        accepted_vars = 0
        discarded_vars = 0

        header = None
        for mu_str in inh:
            if mu_str.startswith('#'):
                continue
            if not header:
                header = mu_str.strip().split()
                continue

            mu_str = mu_str.strip().split()
            variant_annot_types, ensembl_gene_id, mutational_targets = mu_str[0], mu_str[1], mu_str[2]
            mutational_targets = [eval(i) for i in mutational_targets.split(',')]
            variant_annot_types = (variant_annot_types[0], variant_annot_types[1])

            # Filter by specified annotations
            if not variant_annot_types[0] in variant_annots:
                continue
            if suppress_indels_bool and variant_annot_types[1] == 'I':
                continue
            for mu in mutational_targets:
                # Filter by genes accepted in gene_instances
                if not Gene_inst_dict[variant_annot_types].get(ensembl_gene_id):
                    discarded_vars += 1
                    continue
                # gene_instances update
                Gene_inst_dict[variant_annot_types][ensembl_gene_id].mu_list.append(mu)
                accepted_vars += 1

    print(f"Accepted {accepted_vars} variant targets, {discarded_vars} discarded")
    return Gene_inst_dict


def load_muttargs_from_filelist(Gene_inst_dict,
                                input_file_list,
                                variant_annots,
                                suppress_indels_bool):
    """
    Called if --metadata_run_mode is enabled by cohort recurrence methods RaMeDiES-DN or RaMeDiES-CH

    :param Gene_inst_dict: variant annotation type (e.g., CI) -> Ensembl Gene ID -> Gene (object)
    :param input_file_list: comma-separated list of summary statistics input file IDs (initially produced by the
                      "write_muttargs" function), specified by user using the --M parameter
    :param variant_annots: variant annotation types included, specified by user using the --variant_annots parameter
    :param suppress_indels_bool: boolean indicating whether to include indels (False), specified by user
    :return: updated (previously initialized) "Gene_inst_dict" used as an argument to below function "write_muttargs"
    """

    input_file_list = input_file_list.split(',')
    for input_file_id in input_file_list:
        Gene_inst_dict = load_muttargs_from_file(Gene_inst_dict,
                                                 input_file_id,
                                                 variant_annots,
                                                 suppress_indels_bool)
    return Gene_inst_dict


def write_muttargs(Gene_inst_dict, 
                   outfile_mask, 
                   input_directory="<unspecified>", 
                   args_obj=None):
    """
    :param Gene_inst_dict: variant annotation type (e.g., CI) -> Ensembl Gene ID -> Gene (object)
    :param outfile_mask: outfile file prefix, specified by user using the --o parameter
    :param input_directory: directory containing all tab-delimited input variant files per individual
    :param args_obj: dictionary of parsed command-line arguments
    :return: None, write mutational targets by variant type to outfile
    """

    outfile_name = f"{outfile_mask}_{cfg.muttargs_list_DN_ID}.txt"
    with open(outfile_name, 'w') as outh:
        outh.write('# Mutational targets computed from input variant files located in: '+input_directory+'\n')
        init_objs_lib.write_run_info(outh, args_obj)
        outh.write(f"variant_type\tensembl_gene_id\tper_patient_mutational_targets\n")
        for va, va_Gene_obj_dict in Gene_inst_dict.items():
            for ENS_ID, Gene_obj in va_Gene_obj_dict.items():
                muttarg_list = Gene_obj.calculate_muttargs_denovo()
                if len(muttarg_list) == 0:
                    continue

                if args_obj and args_obj.write_muttarg_sums:
                    muttarg_str = np.sum(muttarg_list)
                else:
                    muttarg_str = ','.join([str(i) for i in muttarg_list])
                outh.write(f"{va[0]}{va[1]}\t{ENS_ID}\t{muttarg_str}\n")

    print(f"Variant mutational targets written to: {outfile_name}")


# ==========================================================================================
# Statistics functions used by RaMeDiES-DN (de novo recurrence)
# ==========================================================================================

def count_genes(Gene_inst_dict, return_dict=False):
    """
    :param Gene_inst_dict: variant annotation type (e.g., CI) -> Ensembl Gene ID -> Gene (object)
    :param return_dict: boolean whether to return the complete dictionary (True) or the number of genes (False)
    :return: the number of genes under consideration (or a dictionary of gene_id -> True of the genes under consideration
             consideration for the variant types in question). Used as the Bonferroni correction factor for
             cohort-level recurrence functions RaMeDiES-DN and RaMeDiES-CH.
    """

    ensg_gene_id_dict = {}
    for va in init_objs_lib.variant_type_iter():
        for ensg_gene_id in Gene_inst_dict[va].keys():
            ensg_gene_id_dict[ensg_gene_id] = True
    if return_dict:
        return ensg_gene_id_dict

    return len(ensg_gene_id_dict.keys())


def y_from_mu(gene_mus):
    """
    Called by count_y function

    :param gene_mus: dictionary of Ensembl gene ID -> array of mutational targets (across a cohort)
    :return: sum of observed mutational targets per gene across a cohort
    """

    gene_ys = {}  # Ensembl gene ID -> sum of observed mutational targets across all individuals
    for ensg_gene_id, muttarg_arr in gene_mus.items():
        if muttarg_arr == []:
            continue
        if not gene_ys.get(ensg_gene_id):
            gene_ys[ensg_gene_id] = 0

        gene_ys[ensg_gene_id] += np.sum(muttarg_arr)

    return gene_ys


def count_y(Gene_inst_dict):
    """
    :param Gene_inst_dict: variant annotation type (e.g., CI) -> Ensembl Gene ID -> Gene (object)
    :return: sum of observed mutational targets per gene across a cohort
    """

    gene_mus = {}  # gene_mus: ENSEMBL_ID -> array of mutational tergets
    for va, va_Gene_obj_dict in Gene_inst_dict.items():
        for ensg_gene_id, Gene_obj in va_Gene_obj_dict.items():
            if not gene_mus.get(ensg_gene_id):
                gene_mus[ensg_gene_id] = []
            gene_mus[ensg_gene_id] += Gene_obj.mu_list

    return y_from_mu(gene_mus)  # gene_ys: ENSEMBL ID -> y statistic


def make_gene_mutdict(Gene_inst_dict):
    """
    :param Gene_inst_dict: variant annotation type (e.g., CI) -> Ensembl Gene ID -> Gene (object)
    :return: dictionary of Ensembl gene ID -> list of observed Variant objects (from input tab-delimited variant files)
    """

    gene_mutdict = {}
    for va, va_Gene_obj_dict in Gene_inst_dict.items():
        for ensg_gene_id, Gene_obj in va_Gene_obj_dict.items():
            for Variant_obj in Gene_obj.vars:
                if not gene_mutdict.get(ensg_gene_id):
                    gene_mutdict[ensg_gene_id] = []
                gene_mutdict[ensg_gene_id].append(Variant_obj)
    return gene_mutdict


def calc_denovo_lambda(varcount_dict, Gene_inst_dict, ensg_gene_id):
    """
    :param varcount_dict: dictionary of inheritance type (e.g., M for maternal) -> variant annotation type
                          (e.g., CI for coding indels) -> total observed variants of this type
    :param Gene_inst_dict: variant annotation type (e.g., CI) -> Ensembl Gene ID -> Gene (object)
    :param ensg_gene_id: Ensembl gene ID
    :return: the expected variant count for this gene (i.e., lambda parameter)
    """

    lambda_parameter = 0
    for va in cfg.var_annot_list:
        if Gene_inst_dict[va].get(ensg_gene_id):

            # For each variant annotation type, multiply the variant count by the mutational target.
            lambda_parameter += varcount_dict["DN"][va] * Gene_inst_dict[va][ensg_gene_id].gene_mu
    return lambda_parameter


def false_diag_rate(lambda_parameter, num_samples):
    """
    Computes the binomial upper bound specified by cfg.false_diag_rate
    To obtain an estimate of the fraction of incorrect diagnoses (i.e., number of patients with a variant in a
    recurrently-hit gene, where that variant in that gene in that patient is NOT a true diagnosis), divide the
    false diagnosis rate (computed here as the binomial upper bound specified by cfg.false_diag_rate) by the number of
    variants occurring in a given gene across the cohort.

    :param lambda_parameter: the expected variant count for this gene (i.e., lambda parameter)
    :param num_samples: size of the cohort
    :return: the false diagnosis rate based on the expected mutational target of all variants and the cohort size
    """

    if num_samples == -1:
        return np.nan

    p = lambda_parameter / num_samples  # estimate of the binomial probability

    # true diagnosis rate (i.e., fraction of individuals where this gene is the correct diagnosis)
    higher_p = 1 - cfg.false_diag_rate

    return binom.ppf(higher_p, num_samples, p)


def irwinhall_cdf(x, n):
    """
    :param x: statistic
    :param n: Irwin-Hall parameter (number of elements in a sum of 0-1 Uniforms)
    :return: value of the cumulative density function (CDF) of the Irwin-Hall distribution parameterized by n and
             computed at x
    """

    # Setting CDF = 1 to the right of the realm of definition
    if x >= n:
        return 1

    # Setting CDF = 0 to the left of the realm of definition
    elif x <= 0:
        return 0

    # Central Limit theorem gives a good approximation on about n > 8
    if n > cfg.IH_norm_approx_thr:
        return norm.cdf(x, loc=n / 2, scale=np.sqrt(n / 12))

    k = int(x)
    stat = 0
    for i in range(k + 1):
        xmi = x - i
        if xmi == 0:
            continue

        if i == 0:
            log_factorial_i = 0
        else:
            log_factorial_i = np.sum(np.log([j for j in range(1, i + 1)]))
        if n - i == 0:
            log_factorial_nmi = 0
        else:
            log_factorial_nmi = np.sum(np.log([j for j in range(1, (n - i) + 1)]))

        # Using factorial transform here to make the calculations more robust
        log_prod_xmi = n * np.log(x - i)
        val = ((-1) ** i) * np.exp(-log_factorial_i - log_factorial_nmi + log_prod_xmi)
        stat += val
    return stat


def denovo_prob(lambda_parameter):
    """
    :param lambda_parameter: expected variant count, see "calc_denovo_lambda"
    :return: probability of observing a variant at all given the expected variant count in a cohort
    """

    return 1.0 - np.exp(-lambda_parameter)


def process_single_gene(y, lambda_parameter, gene_constraint_weight, num_samples):
    """
    For more information on the weighted FDR procedure, see Genovese et al., 2006 for details
    (https://www.jstor.org/stable/20441304)

    :param y: sum of mutational targets of independent variants within a single gene
    :param lambda_parameter: (float) expected variant count in this gene
    :param gene_constraint_weight: gene constraint based weight with mean of 1 used in the weighted FDR procedure
    :param num_samples: number of samples in the cohort
    :return: the Irwin-Hall/Poisson recurrence statistic computed for a single gene
    """

    P_dnv = denovo_prob(lambda_parameter)
    P_val = 0

    # The "infinite" sum is computed with the first ~1000 values
    for i in range(1, cfg.maxIHval):

        # Irwin-Hall survival function
        P_IH = 1 - irwinhall_cdf(y, i)

        # Very small values may be negative due to computational errors
        if P_IH < 0:
            P_IH = 0

        P_Pois = poisson.pmf(i, lambda_parameter)

        # Element of the sum
        add_P = P_IH * P_Pois

        # Early stop in case of the added value being small with respect to the already calculated value
        if P_val != 0 and add_P / P_val < cfg.pval_precision:
            break

        P_val += add_P

    # Probability of an observed mutational pattern given that a single mutation has been observed. Used in Q-Q plots
    P_cond = P_val / P_dnv

    # Q value of the weighted FDR procedure
    if gene_constraint_weight != 0:
        Q_val = P_val / gene_constraint_weight
    else:
        Q_val = 1.0

    # Computing the false diagnosis rate
    f_diag_rate = false_diag_rate(lambda_parameter, num_samples)

    return [Q_val, P_val, P_cond, P_dnv, lambda_parameter, f_diag_rate, gene_constraint_weight]


def generate_output_line(pvalues_array,
                         out_handle,
                         ensembl_gene_id,
                         ensembl_to_genename,
                         gene_mutdict,
                         first_line=False):
    """
    :param pvalues_array: output of "process_single_gene", array of [Q_val, P_val, P_cond, P_dnv, lambda_parameter,
                          false_diagnosis_rate, gene_constraint_weight]
    :param out_handle: output file object (opened for write)
    :param ensembl_gene_id: Ensembl gene ID
    :param ensembl_to_genename: dictionary of ensembl gene ID -> HGNC gene name
                                (see init_objs_lib/make_ENS2GeneID_dict for details) for details)
    :param gene_mutdict: output of "make_gene_mutdict", array of Variant objects corresponding to variants present in
                         tab-delimited input variant files
    :param first_line: boolean indicating whether to write a header (True)
    :return: None, but write properly formatted lines to open file handle as specified
    """

    if first_line:
        out_handle.write('# Genes harboring de novo mutations across cohort ranked by Q-value\n')
        out_handle.write('# file_names should include relevant patient identifiers\n')
        out_handle.write('# Variant inheritance is M=maternally-inherited, P=paternally-inherited, and DN=denovo\n')
        out_handle.write('# variant_info is &-delimited values: chromsome:refallele_position_altallele|variant_type|' +
                         'variant_functionality_score|variant_inheritance|file_name\n')

        print_arr = ["file_names",
                     "ensembl_gene_id",
                     "gene_name",
                     "Q_val",
                     "P_val",
                     "P_cond",
                     "P_dnv",
                     "poisson_lambda",
                     "false_diag_rate",
                     "gene_weight",
                     "variant_info"]

        out_handle.write('\t'.join(print_arr) + '\n')
        return None

    # DEFAULT run: variant and individual information is included in output
    if gene_mutdict:

        # List of individuals carrying detected mutations in the gene
        ind_list = [Variant_obj.proband_id for Variant_obj in gene_mutdict[ensembl_gene_id]]
        ind_list = ','.join(ind_list)

        # List of variant information
        var_info_list = [Variant_obj.print_info() for Variant_obj in gene_mutdict[ensembl_gene_id]]
        var_info_list = ','.join(var_info_list)

    # METADATA run: variant- and individual-level data is NOT printed
    else:
        ind_list = '.'
        var_info_list = '.'

    gene_name = ensembl_to_genename.get(ensembl_gene_id)

    if not gene_name:
        gene_name = '.'

    print_arr = [ind_list,
                 ensembl_gene_id,
                 gene_name,
                 str(pvalues_array[0]),
                 str(pvalues_array[1]),
                 str(pvalues_array[2]),
                 str(pvalues_array[3]),
                 str(pvalues_array[4]),
                 str(pvalues_array[5]),
                 str(pvalues_array[6]),
                 var_info_list]

    out_handle.write('\t'.join(print_arr) + '\n')


def weighted_fdr_correction(qval_array, fdrs_to_report, num_genes):
    """
    Perform Benjamini-Hochberg algorithm on Q-values for weighted FDR procedure
    For details, see Genovese et al., 2006 (https://www.jstor.org/stable/20441304)

    :param qval_array: array of unique Q-values
    :param fdrs_to_report: FDR limits at which to report output, specified by user, default 5% and 10%
    :param num_genes: number of genes processed
    :return: the Q-value threshold corresponding to each FDR threshold specified in "fdrs_to_report"
    """

    qval_array = sorted(qval_array)
    qval_array_length = len(qval_array)

    fdr_cutoff_to_qval_threshold = {fdr: 0 for fdr in fdrs_to_report}

    # For each FDR cutoff to report, compute the corresponding Q-value threshold
    for i in range(qval_array_length):
        for fdr in fdrs_to_report:
            qval_threshold = fdr * (i + 1) / num_genes
            if qval_array[i] < qval_threshold:
                fdr_cutoff_to_qval_threshold[fdr] = qval_threshold

    return fdr_cutoff_to_qval_threshold


def FDR_procedure(outfile_mask, fdrs_to_report, num_genes):
    """
    :param outfile_mask: outfile file prefix, specified by user using the --o parameter
    :param fdrs_to_report: FDR limits at which to report output, specified by user, default 5% and 10%
    :param num_genes: number of genes processed
    :return: None, but write out intermediate file used to perform the weighted FDR procedure on P-values
    """

    outfile_dict = {}  # qvalue -> output string; used to produce output sorted by qvalue

    qval_array = []

    filename = f"{outfile_mask}_{cfg.DN_result}_tmp.txt"

    with open(filename, 'r') as inh:
        header = None
        for output_str in inh:
            if output_str.startswith('#'):
                continue
            if not header:
                header = {k: i for i, k in enumerate(output_str.strip().split())}
                header_line = output_str.strip().split()
                continue

            output_str = output_str.strip().split('\t')
            qval = eval(output_str[header["Q_val"]])
            qval_array.append(qval)
            if not outfile_dict.get(qval):
                outfile_dict[qval] = []
            outfile_dict[qval].append(output_str)

    fdr_to_qval_threshold = weighted_fdr_correction(qval_array, fdrs_to_report, num_genes)
    fdr_column_names = '\t'.join([f"FDR_{p}" for p in fdr_to_qval_threshold])

    # Final output file
    filename = f"{outfile_mask}_{cfg.DN_result}.txt"
    with open(filename, 'w') as outh:

        # write out original header plus additional FDR-relevant columns
        hs = '\t'.join(header_line)
        outh.write(f"{hs}\t{fdr_column_names}\n")

        for qval in sorted(outfile_dict.keys()):
            fdr_results_array = []

            # Using Q-value thresholds from weighted_fdr_correction
            for FDR in fdr_to_qval_threshold:

                # FDR pass/ not pass is added to the output
                if qval < fdr_to_qval_threshold[FDR]:
                    fdr_results_array.append("True")
                else:
                    fdr_results_array.append("False")

            concatenated_fdr_results = '\t'.join(fdr_results_array)
            for file_str_arr in outfile_dict[qval]:
                fsa = '\t'.join(file_str_arr)
                outh.write(f"{fsa}\t{concatenated_fdr_results}\n")

    print(f"Results written to {filename}")


def calc_denovo_stat(ensg_to_y_stat,
                     varcount_dict,
                     Gene_inst_dict,
                     outfile_mask,
                     num_samples,
                     gene_mutdict,
                     fdrs_to_report_str,
                     gene_constraint_score_name,
                     num_genes):
    """
    Master function for RaMeDiES-DN; based on input variant information OR supplied metadata, calculates
    a de novo recurrence p-value per gene, performs a weighted FDR procedure using gene constraint weights,
    and writes results to file.

    :param ensg_to_y_stat: dictionary of Ensembl gene ID -> (y statistic, count_y output)
    :param varcount_dict: inheritance type (e.g., M = maternal) -> variant annotation type (e.g., CI for coding indels)
                          -> # variants of this type
    :param Gene_inst_dict: variant annotation type (e.g., CI) -> Ensembl Gene ID -> Gene (object)
    :param outfile_mask: outfile file prefix, specified by user using the --o parameter
    :param num_samples: number of samples in cohort (with corresponding input variant files)
    :param gene_mutdict: output of "make_gene_mutdict"; Ensembl gene ID -> list of Variant objects extracted from
                         input tab-delimited variant files
    :param fdrs_to_report_str: comma-separated list of FDR thresholds (specified by user)
    :param gene_constraint_score_name: name of the gene-based constraint score to use (e.g., GeneBayes)
    :param num_genes: number of processed genes
    :return: None, but perform all operations and create output files as required by RaMeDiES-DN
    """

    # FDR values are obtained from the input string
    fdrs_to_report = [eval(i) for i in fdrs_to_report_str.split(',')]

    # constraint_bins: Gene ENSEMBL ID -> gene constraint enrichment bin
    # Enrichment values are the proportions of exclusively dominant OMIM genes in gene constraint score deciles
    constraint_bins = init_objs_lib.load_gene_score_bins(gene_constraint_score_name)

    # ensg_to_genename: ensembl_gene_id -> Gene ID
    ensg_to_genename = init_objs_lib.make_ENS2GeneID_dict()

    # Intermediate output file
    filename = f"{outfile_mask}_{cfg.DN_result}_tmp.txt"

    with open(filename, 'w') as out_handle:
        # Header is printed
        generate_output_line(pvalues_array=None,
                             out_handle=out_handle,
                             ensembl_gene_id=None,
                             ensembl_to_genename=ensg_to_genename,
                             gene_mutdict=gene_mutdict,
                             first_line=True)

        # Iterating over genes
        for ensembl_gene_id, y in ensg_to_y_stat.items():

            # Expected number of variant calculated
            lambda_parameter = calc_denovo_lambda(varcount_dict, Gene_inst_dict, ensembl_gene_id)

            # s_het weight is obtained
            gene_constraint_weight = constraint_bins.get(ensembl_gene_id)

            # For gene with no s_het values, s_het weight is set to unity
            if not gene_constraint_weight:
                gene_constraint_weight = 1.0

            # Calculating the probabilities and Q values
            # See process_single_gene for details
            P_arr = process_single_gene(y, lambda_parameter, gene_constraint_weight, num_samples)

            # Intermediate output is printed
            generate_output_line(pvalues_array=P_arr,
                                 out_handle=out_handle,
                                 ensembl_gene_id=ensembl_gene_id,
                                 ensembl_to_genename=ensg_to_genename,
                                 gene_mutdict=gene_mutdict,
                                 first_line=False)

    # Weighted FDR step
    FDR_procedure(outfile_mask, fdrs_to_report, num_genes)

    # Remove the intermediate output
    remove(filename)


# ==========================================================================================
# Write and load summary statistics files containing comphet VARIANT COUNT information
# ==========================================================================================

def write_mutnum_prods(mutnum_prod_dict,
                       mutnum_prod_dist_dict,
                       outfile_mask,
                       args_obj,
                       input_directory="<unspecified>"):
    """
    Called by "mutnum_prod" function below. Write mutation number products and distribution to file

    :param mutnum_prod_dict: pair of variant annotations (specified in cfg.var_annot_list) ->
                             sum of the products of their counts. See "mutnum_prod" function below for details.
    :param mutnum_prod_dist_dict: pair of variant annotations -> product value -> frequency. Used for QC purposes.
    :param outfile_mask: outfile file prefix, specified by user using the --o parameter
    :param args_obj: dictionary of parsed command-line arguments
    :param input_directory: directory containing all tab-delimited input variant files per individual
    :return: None, but write the mutation number products to file
    """

    filename = f"{outfile_mask}_{cfg.mutnum_prod_CH}.txt"

    # Writing the mutation products to the intermediate output
    with open(filename, 'w') as outh:
        outh.write("# Sum (across cohort) of all products of (# paternally-inherited variants)x(# maternally-inherited variants) for each variant type pair\n")
        outh.write('# Variant input files processed from: '+input_directory+'\n')
        init_objs_lib.write_run_info(outh, args_obj)
        outh.write("paternal_variant_type\tmaternal_variant_type\ttotal_product_of_variant_counts\n")
        for annot_pair, count in mutnum_prod_dict.items():
            annot_P = f"{annot_pair[0][0]}{annot_pair[0][1]}"
            annot_M = f"{annot_pair[1][0]}{annot_pair[1][1]}"
            outh.write(f"{annot_P}\t{annot_M}\t{count}\n")

    print(f"Variant count products written to {filename}")

    filename = f"{outfile_mask}_{cfg.mutnum_prod_dist_CH}.txt"

    # Writing the distribution of mutation products to the intermediate output
    with open(filename, 'w') as outh:
        outh.write('# Number of input variant files with specific (paternal_variant_count)x(maternal_variant_count) values\n')
        outh.write('# Variant input files processed from: '+input_directory+'\n')
        init_objs_lib.write_run_info(outh, args_obj)
        outh.write("paternal_variant_type\tmaternal_variant_type\tproduct_of_variant_counts\tnumber_samples\n")
        for annot_pair, dist_dict in mutnum_prod_dist_dict.items():
            annot_P = f"{annot_pair[0][0]}{annot_pair[0][1]}"
            annot_M = f"{annot_pair[1][0]}{annot_pair[1][1]}"
            for mnp in sorted(dist_dict.keys()):
                freq = dist_dict[mnp]
                outh.write(f"{annot_P}\t{annot_M}\t{mnp}\t{freq}\n")

    print(f"Distributions of variant count products written to {filename}")


def mutnum_prod(varcount_dict, outfile_mask, args_obj, input_directory="<unspecified>"):
    """
    Products of mutation counts are used to compute the expected numbers of compound heterozygous variants,
    i.e., the Poisson lambda parameter

    :param varcount_dict: proband_id (filename) -> inheritance type -> variant annotation type -> number of variants
    :param outfile_mask: outfile file prefix, specified by user using the --o parameter
    :param args_obj: dictionary of parsed command-line arguments
    :param input_directory: directory containing all tab-delimited input variant files per individual
    :return: dictionary pair of variant annotations (specified in cfg.var_annot_list) ->
             sum of the products of their counts. See "mutnum_prod" function below for details.
    """

    mutnum_prod_dict = {}  # pair of variant annotations -> sum of the products of their counts
    mutnum_prod_dist_dict = {}  # pair of variant annotations -> products of their counts -> frequency

    # Initiating the double variant count dictionaries
    # cfg.var_annot_list: list of variant annotations
    # 	1st letter: 'C' for coding, 'I' for intronic
    # 	2nd letter: 'S' for SNP, 'I' for indel
    for va1 in cfg.var_annot_list:
        for va2 in cfg.var_annot_list:
            mutnum_prod_dict[(va1, va2)] = 0
            mutnum_prod_dist_dict[(va1, va2)] = {}

    # Iterating over individuals, summing the mutation count products
    for proband_id, inher_dict in varcount_dict.items():
        if not inher_dict.get('M') or not inher_dict.get('P'):
            continue

        # Iterating over the pairs of annotations
        for annot_pair in mutnum_prod_dict.keys():
            va_P = annot_pair[0]
            va_M = annot_pair[1]
            mnp = inher_dict['P'][va_P] * inher_dict['M'][va_M]
            mutnum_prod_dict[(va_P, va_M)] += mnp
            if not mutnum_prod_dist_dict[(va_P, va_M)].get(mnp):
                mutnum_prod_dist_dict[(va_P, va_M)][mnp] = 0
            mutnum_prod_dist_dict[(va_P, va_M)][mnp] += 1

    # Writing intermediate outputs, which are used in the metadata runs
    write_mutnum_prods(mutnum_prod_dict,
                       mutnum_prod_dist_dict,
                       outfile_mask,
                       args_obj,
                       input_directory)

    return mutnum_prod_dict


def read_mutnum_prods_from_file(file_id,
                                mutnum_prod_dict,
                                variant_annots,
                                suppress_indels_bool):
    """
    Called by "read_mutnum_prods" and processes a single file in the case of a metadata run on summary statistics

    :param file_id: summary statistics input file ID (initially produced by the "write_mutnum_prods" function),
                     specified by user using the --M parameter
    :param mutnum_prod_dict: dictionary of pair of variant annotations -> sum of the products of their counts
    :param variant_annots: variant annotation types included, specified by user using the --variant_annots parameter
    :param suppress_indels_bool: boolean indicating whether to include indels (False), specified by user
    :return: updated dictionary of pair of variant annotations -> sum of the products of their counts
    """

    filename = f"{file_id}_{cfg.mutnum_prod_CH}.txt"
    with open(filename, 'r') as inh:
        for mnp_str in inh:

            header = None
            if mnp_str.startswith('#'):
                continue
            if not header:
                header = mnp_str.strip().split()
                continue

            mnp_str = mnp_str.strip().split()

            va_P = (mnp_str[0][0], mnp_str[0][1])  # Paternal variant annotation

            va_M = (mnp_str[1][0], mnp_str[1][1])  # Maternal variant annotation

            # Check if the line meets run specifications
            if va_P[0] not in variant_annots or va_M[0] not in variant_annots:
                continue

            # Specified regime of indel processing
            if suppress_indels_bool and (va_P[1] == 'I' or va_M[1] == 'I'):
                continue

            count = eval(mnp_str[2])
            mutnum_prod_dict[(va_P, va_M)] += count

    return mutnum_prod_dict


def read_mutnum_prods(file_id_list, variant_annots, suppress_indels_bool):
    """
    :param file_id_list: comma-separated list of summary statistics input file IDs (initially produced by the
                         "write_mutnum_prods" function), specified by user using the --M parameter
    :param variant_annots: variant annotation types included, specified by user using the --variant_annots parameter
    :param suppress_indels_bool: boolean indicating whether to include indels (False), specified by user
    :return: dictionary of pair of variant annotations -> sum of the products of their counts
    """

    mutnum_prod_dict = {}  # pair of variant annotations -> sum of the products of their counts

    # Initiating the double variant count dictionaries
    # cfg.var_annot_list: list of variant annotations
    # 	1st letter: 'C' for coding, 'I' for intronic
    # 	2nd letter: 'S' for SNP, 'I' for indel
    for va1 in cfg.var_annot_list:
        for va2 in cfg.var_annot_list:
            mutnum_prod_dict[(va1, va2)] = 0

    for file_id in file_id_list.split(','):

        # For each specified ID, call read_mutnum_prods_from_file
        mutnum_prod_dict = read_mutnum_prods_from_file(file_id,
                                                       mutnum_prod_dict,
                                                       variant_annots,
                                                       suppress_indels_bool)

    return mutnum_prod_dict


# ==========================================================================================
# Write and load summary statistics files containing comphet MUTATIONAL TARGET information
# ==========================================================================================

def write_comphet_muttargs(gene_comphet_mutdict, outfile_mask, input_directory="<unspecified>"):
    """
    :param gene_comphet_mutdict: Ensembl gene ID -> list of CH_Variant objects extracted from input tab-delimited 
                                 variant files
    :param outfile_mask: outfile file prefix, specified by user using the --o parameter
    :param input_directory: directory containing all tab-delimited input variant files per individual
    :return: None, but write out mutational target information to file
    """

    filename = f"{outfile_mask}_{cfg.muttargs_CH_ID}.txt"
    with open(filename, 'w') as outh:
        outh.write(
            '# Compound heterozygous mutational targets computed for processed variant files in ' + input_directory + '\n')
        outh.write('# variant types are in the format paternal_variant_type,maternal_variant_type\n')
        outh.write(f"variant_type\tensembl_gene_id\tper_patient_mutational_targets\n")
        for ensembl_gene_id, VariantCollection_obj in gene_comphet_mutdict.items():
            for CH_variant_obj in VariantCollection_obj.CH_list:

                # va_P: variant annotation of the paternally inherited variant
                # va_M: variant annotation of the maternally inherited variant
                va_P = CH_variant_obj.var_P.var_annot
                va_M = CH_variant_obj.var_M.var_annot

                # CH annotation
                va = f"{va_P[0]}{va_P[1]},{va_M[0]}{va_M[1]}"

                # CH mutational target
                muttarg = round(CH_variant_obj.mu2, 3)
                outh.write(f"{va}\t{ensembl_gene_id}\t{muttarg}\n")

    print(f"Compound heterozygous variant mutational targets written to: {filename}")


def read_comphet_muttargs_from_file(file_id,
                                    CH_y_dict,
                                    variant_annots,
                                    suppress_indels_bool):
    """
    Called by "read_comphet_muttargs"

    :param file_id: summary statistics input file ID (initially produced by the "write_mutnum_prods" function),
                    specified by user using the --M parameter
    :param CH_y_dict: ensembl gene ID -> comphet y statistic (initialized in read_comphet_muttargs, updated here)
    :param variant_annots: variant annotation types included, specified by user using the --variant_annots parameter
    :param suppress_indels_bool: boolean indicating whether to include indels (False), specified by user
    :return: updated ensembl gene ID -> comphet y statistic
    """

    filename = f"{file_id}_{cfg.muttargs_CH_ID}.txt"
    print("Reading compound heterozygous mutational targets from file:", filename)
    with open(filename, 'r') as inh:
        header = None
        for mt_line in inh:
            if mt_line.startswith('#'):
                continue
            if not header:
                header = mt_line.strip().split()
                continue

            mt_line = mt_line.strip().split()
            variant_type = mt_line[0].split(',')
            ensembl_gene_id = mt_line[1]
            mutational_targets = [eval(i) for i in mt_line[2].split(',')]

            # only consider the variant types (coding/intronic) that are specified:
            if variant_type[0][0] not in variant_annots or variant_type[1][0] not in variant_annots:
                continue
            if suppress_indels_bool and (variant_type[0][1] == 'I' or variant_type[1][1] == 'I'):
                continue

            if not CH_y_dict.get(ensembl_gene_id):
                CH_y_dict[ensembl_gene_id] = 0

            for muttarg in mutational_targets:
                CH_y_dict[ensembl_gene_id] += 1 - muttarg

    return CH_y_dict


def read_comphet_muttargs(file_id_list, variant_annots, suppress_indels_bool):
    """
    In the case of a metadata run, read the comphet mutational target information from the list of files
    originally produced by "write_comphet_muttargs"

    :param file_id_list: comma-separated list of summary statistics input file IDs (initially produced by the
                         "write_mutnum_prods" function), specified by user using the --M parameter
    :param variant_annots: variant annotation types included, specified by user using the --variant_annots parameter
    :param suppress_indels_bool: boolean indicating whether to include indels (False), specified by user
    :return: ensembl gene ID -> comphet y statistic
    """

    CH_y_dict = {}

    for file_id in file_id_list.split(','):
        CH_y_dict = read_comphet_muttargs_from_file(file_id,
                                                    CH_y_dict,
                                                    variant_annots,
                                                    suppress_indels_bool)

    return CH_y_dict


# ==================================================================================
# Statistics functions used by RaMeDiES-CH (compound heterozygous recurrence)
# ==================================================================================

def calc_comphet_lambda(mutnum_prod_dict, Gene_inst_dict, ensembl_gene_id):
    """
    :param mutnum_prod_dict: dictionary of pair of variant annotations -> sum of the products of their counts
    :param Gene_inst_dict: variant annotation type (e.g., CI) -> Ensembl Gene ID -> Gene (object)
    :param ensembl_gene_id: Ensembl gene ID
    :return: expected number of comphet variants in a gene given the cohort size; Poisson lambda parameter
    """

    comphet_lambda = 0

    for (va_P, va_M), count in mutnum_prod_dict.items():
        if not Gene_inst_dict[va_P].get(ensembl_gene_id) or not Gene_inst_dict[va_M].get(ensembl_gene_id):
            continue
        mu_P = Gene_inst_dict[va_P][ensembl_gene_id].gene_mu
        mu_M = Gene_inst_dict[va_M][ensembl_gene_id].gene_mu

        # The expected number is updated based on the gene mutational targets
        #	and the products of mutation counts
        comphet_lambda += count * mu_P * mu_M

    return comphet_lambda


def calc_CH_y_dict(gene_comphet_mutdict):
    """
    :param gene_comphet_mutdict: Ensembl gene ID -> list of CH_Variant objects extracted from input tab-delimited
                                 variant files
    :return: updated dictionary of ensembl gene ID -> comphet y statistic
    """

    CH_y_dict = {}
    for ensembl_gene_id, VariantCollection_obj in gene_comphet_mutdict.items():
        if VariantCollection_obj.CH_list == []:
            continue

        y = VariantCollection_obj.calc_comphet_y()
        CH_y_dict[ensembl_gene_id] = y

    return CH_y_dict


def generate_output_line_comphet(pvalues_array,
                                 out_handle,
                                 gene_comphet_mutdict,
                                 ensembl_gene_id,
                                 ensg_to_genename,
                                 first_line):
    """
    :param pvalues_array: output of "process_single_gene", array of [P_val, P_cond, P_comphet, lambda_parameter,
                          false_diagnosis_rate]
    :param out_handle: output file object (opened for write)
    :param gene_comphet_mutdict: output of "make_gene_comphet_mutdict", array of Variant objects corresponding to
                                 variants present in tab-delimited input variant files
    :param ensembl_gene_id: Ensembl gene ID
    :param ensg_to_genename: dictionary of ensembl gene ID -> HGNC gene name
                                (see init_objs_lib/make_ENS2GeneID_dict for details) for details)
    :param first_line: boolean indicating whether to write a header (True)
    :return: None, but write properly formatted lines to open file handle as specified
    """

    if first_line:
        # Prints the first line
        print_arr = ["file_names",
                     "ensembl_gene_id",
                     "gene_name",
                     "P_val",
                     "P_cond",
                     "P_comphet",
                     "poisson_lambda",
                     "false_diagnosis_rate",
                     "variant_info"]

        out_handle.write('\t'.join(print_arr) + '\n')
        return None

    # DEFAULT run: individual variant-level information is printed to output
    if gene_comphet_mutdict:

        # IDs of individuals are printed
        ind_list = [CH.var_P.proband_id for CH in gene_comphet_mutdict[ensembl_gene_id].CH_list]
        ind_list = ','.join(ind_list)

        # Variant information is printed
        var_info_list = [CH.print_info() for CH in gene_comphet_mutdict[ensembl_gene_id].CH_list]
        var_info_list = ','.join(var_info_list)

    # METADATA run: all individual-level information is excluded from output file
    else:
        ind_list = "."
        var_info_list = '.'

    gene_name = ensg_to_genename[ensembl_gene_id]

    print_arr = [ind_list,
                 ensembl_gene_id,
                 gene_name,
                 str(pvalues_array[0]),
                 str(pvalues_array[1]),
                 str(pvalues_array[2]),
                 str(pvalues_array[3]),
                 str(pvalues_array[4]),
                 var_info_list]

    return '\t'.join(print_arr) + '\n'


def calc_comphet_stat(CH_y_dict,
                      gene_comphet_mutdict,
                      Gene_inst_dict,
                      mutnum_prod_dict,
                      outfile_mask,
                      num_samples,
                      num_genes):
    """
    Master function for RaMeDiES-CH; based on input variant information OR supplied metadata, calculates
    a comphet recurrence p-value per gene and writes results to file.

    :param CH_y_dict: ensembl gene ID -> comphet y statistic
    :param gene_comphet_mutdict: output of "make_gene_comphet_mutdict", array of Variant objects corresponding to
                                 variants present in tab-delimited input variant files
    :param Gene_inst_dict: variant annotation type (e.g., CI) -> Ensembl Gene ID -> Gene (object)
    :param mutnum_prod_dict: dictionary of pair of variant annotations -> sum of the products of their counts
    :param outfile_mask: outfile file prefix, specified by user using the --o parameter
    :param num_samples: number of samples in cohort
    :param num_genes: number of genes under consideration
    :return: None, but write properly formatted output lines as specified
    """

    ensg_to_genename = init_objs_lib.make_ENS2GeneID_dict()

    # Output file
    filename = f"{outfile_mask}_{cfg.CH_result}.txt"

    # P-value -> list of output lines. Used to make a sorted output file
    output_dict = {}

    with open(filename, 'w') as out_handle:
        # Printing the output header
        generate_output_line_comphet(pvalues_array=None,
                                     out_handle=out_handle,
                                     gene_comphet_mutdict=gene_comphet_mutdict,
                                     ensembl_gene_id=None,
                                     ensg_to_genename=ensg_to_genename,
                                     first_line=True)

        # For each gene, calculate statistic and put the output line into output_dict
        for ensembl_gene_id, y in CH_y_dict.items():
            comphet_lambda_parameter = calc_comphet_lambda(mutnum_prod_dict, Gene_inst_dict, ensembl_gene_id)
            pvalues_array = process_single_gene(y,
                                        comphet_lambda_parameter,
                                        gene_constraint_weight=1.0,  # Dummy value
                                        num_samples=num_samples)[1:-1]  # Discarding Q and gene constraint values

            ol = generate_output_line_comphet(pvalues_array,
                                              out_handle,
                                              gene_comphet_mutdict,
                                              ensembl_gene_id,
                                              ensg_to_genename,
                                              first_line=False)

            if not output_dict.get(pvalues_array[0]):
                output_dict[pvalues_array[0]] = []
            output_dict[pvalues_array[0]].append(ol)

        # Writing the otput data
        for P in sorted(output_dict.keys()):
            for ol in output_dict[P]:
                out_handle.write(ol)

    print(f"Results written to {filename}")


# ==================================================================================
# Statistics functions used by RaMeDiES-IND (individual-level compound heterozygous variants)
# ==================================================================================

def clean_Gene_inst_dict(Gene_inst_dict):
    """
    Clear out variant information from existing gene instance dictionary to reuse the initially-built
    dictionary multiple times

    :param Gene_inst_dict: variant annotation type (e.g., CI) -> Ensembl Gene ID -> Gene (object)
    :return: a new Gene_inst_dict dictionary with all variants and mutational targets removed
    """
    # gene_instances: variant annotation -> ensembl_gene_id -> Gene object
    for va in Gene_inst_dict.keys():
        for ensembl_gene_id in Gene_inst_dict[va].keys():
            # vars and mu_list attributes are cleaned
            Gene_inst_dict[va][ensembl_gene_id].vars = []
            Gene_inst_dict[va][ensembl_gene_id].mu_list = []
    return Gene_inst_dict


def calc_comphet_individual_pvalue(y, comphet_lambda_parameter):
    """
    :param y: per-individual comphet test statistic, which is the fraction of comphet mutational targets across the
              genome that is equal to or lower than the observed mutational target
    :param comphet_lambda_parameter: expected number of comphets in this individual given observed variant counts
    :return: array of [gene_pvalue, conditional_pvalue (for QQ plots), same starting y statistic, same starting
                       lambda parameter]
    """

    gene_pvalue = 0
    stat_val = 1
    K_i = 1

    # Infinite sum approximation
    while stat_val > 1e-16 or K_i < 5 * comphet_lambda_parameter or K_i < 10:
        stat_val = (1 - (1 - y) ** K_i) * poisson.pmf(K_i, comphet_lambda_parameter)
        gene_pvalue += stat_val
        K_i += 1

    # Probability of y conditional on observation of at least one CH
    conditional_pvalue = gene_pvalue / (1 - np.exp(-comphet_lambda_parameter))
    return [gene_pvalue, conditional_pvalue, y, comphet_lambda_parameter]


def print_comphet_individual_line(out_handle,
                                  pvalues_array=None,
                                  top_comphet=None,
                                  ensembl_gene_id=None,
                                  ensg_to_genename=None,
                                  first_line=True):
    """
    :param out_handle: output file object (opened for write)
    :param pvalues_array: output of "calc_comphet_individual_pvalue", array of [P_val, P_cond, P_comphet, lambda_parameter,
                          false_diagnosis_rate]
    :param top_comphet: CH_variant object corresponding to the least expected comphet occurrence
    :param ensembl_gene_id: Ensembl gene ID
    :param ensg_to_genename: dictionary of ensembl_gene_id -> HGNC gene name
    :param first_line: boolean indicating whether only the header should be written (True) or only value lines (False)
    :return: None, but write properly formatted lines to open file handle as specified
    """

    if first_line:
        print_arr = ["file_name",
                     "ensembl_gene_id",
                     "gene_name",
                     "P_val",
                     "P_cond",
                     "y_stat",
                     "poisson_lambda",
                     "variant_info"]

        out_handle.write('\t'.join(print_arr) + '\n')
        return None

    # Print the output for a single individual
    gene_name = ensg_to_genename[ensembl_gene_id]
    var_info = top_comphet.print_info()
    print_arr = [top_comphet.var_P.proband_id,
                 ensembl_gene_id,
                 gene_name,
                 str(pvalues_array[0]),
                 str(pvalues_array[1]),
                 str(pvalues_array[2]),
                 str(pvalues_array[3]),
                 var_info]

    out_handle.write('\t'.join(print_arr) + '\n')


def calc_comphet_individual_statistic(out_handle,
                                      variant_counts,
                                      top_comphet,
                                      Gene_inst_dict,
                                      ensg_to_genename):
    """
    Master function for RaMeDiES-IND; based on input variant information OR supplied metadata, calculates
    a p-value for the most surprising comphet per individual in a cohort and writes to file.

    :param out_handle: output file object (opened for write)
    :param variant_counts: dictionary of inheritance type (e.g., M for maternal) -> variant annotation type
                           (e.g., CI for coding indel) -> # of variants of this type per individual
    :param top_comphet: CH_variant object corresponding to the least expected comphet occurrence
    :param Gene_inst_dict: variant annotation type (e.g., CI) -> Ensembl Gene ID -> Gene (object)
    :param ensg_to_genename: dictionary of ensembl_gene_id -> HGNC gene name
    :return: None, but write properly formatted output lines as specified
    """

    comphet_lambda_parameter = 0  # expected number of comphet variants per individual

    # per-individual comphet test statistic, which is the fraction of comphet mutational targets across the genome
    # that are equal to or lower than the observed comphet mutational target.
    y = 0

    # corr_factor: maximal CH mutational target in the sense of y
    corr_factor = 0  # maximal comphet mutational target (as ranked with "y")

    # Iterate over pairs of variant annotations
    for va_P, var_num_P in variant_counts['P'].items():
        for va_M, var_num_M in variant_counts['M'].items():
            if var_num_P * var_num_M == 0:
                continue

            # Iterate over genes
            for ensembl_gene_id, Gene_obj_P in Gene_inst_dict[va_P].items():
                Gene_obj_M = Gene_inst_dict[va_M].get(ensembl_gene_id)
                if not Gene_obj_M:
                    continue

                # Counting all mutational targets lesser than top_comphet.mu ** 2
                mut_targ_P = min(Gene_obj_P.gene_mu, top_comphet.mu)
                mut_targ_M = min(Gene_obj_M.gene_mu, top_comphet.mu)

                # y, comphet_lambda_parameter and corr_factor are updated here
                # NOTE that the observed mutational target is actually top_comphet.mu ** 2, which might be
                #      greater than the total squared mutational target of the gene, which is why we
                y += np.max([mut_targ_P, mut_targ_M]) ** 2
                corr_factor += np.max([Gene_obj_P.gene_mu, Gene_obj_M.gene_mu]) ** 2
                comphet_lambda_parameter += (var_num_P * Gene_obj_P.gene_mu) * (var_num_M * Gene_obj_M.gene_mu)

    # y normalization
    y = y / corr_factor

    # probability of y given comphet_lambda_parameter
    pvalues_array = calc_comphet_individual_pvalue(y, comphet_lambda_parameter)

    # Writing the output
    print_comphet_individual_line(out_handle,
                                  pvalues_array,
                                  top_comphet,
                                  top_comphet.ENS_ID,
                                  ensg_to_genename,
                                  first_line=False)

# Written on 01.23.2024 by Mikhail Moldovan, HMS DBMI
