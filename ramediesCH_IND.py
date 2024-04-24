"""
ramediesCH_IND

Main script for the RaMeDiES individual-level compound heterozygote analysis

Python version: 3.6 and above

Module requirements:
1. numpy
2. argparse
3. os

Next goes the description of the program. For more information about the usage,
consult the manual on GitHub and the help menu via the command:
python ramediesDN.py -h

######
First, the collection of Gene objects is calculated.
Gene objects are specified in stat_lib.

Collection is stratified by annotations (coding/intronic variants and
SNPs/indels) and by gene ENSEMBL IDs and is organized as a two-key dictionary:

gene_instances: variant annotation -> ENSEMBL ID -> Gene object.

gene_instances is initialized by the parse_variant_scores_files function
from the init_objs_lib.py library. Gene class is specified in the stat_lib
library.

Instances of the Gene class contain arrays of score-mutational targets 
distributions and the information about specific variants of a given annotation 
that landed in a given gene. Variant-level information is specified in the 
Variant class in stat_lib.

######
Second, for each individual, the VariantCollection objects (specified in stat_lib) 
as initiated for each gene. Unlike Gene objects, VariantCollection objects are not 
tied to specific variant annotations, which allows to infer all CH variants by
using the make_CH_var_list method of VariantCollection class and to calculate 
the test statistic (called y) by the calc_CH_y method of VariantCollection.
The products of per-cohort variant counts are also calculated. Finally,
a statistic described in our paper is calculated and the result is written to the
output file with identified provided through the --o option.

If you have any issues, contact us at mikhail_moldovan@hms.harvard.edu
"""

# Libraries 
import init_objs_lib as iol
import cfg
import parse_VCF_lib as pvl
import stat_lib as sl

# Standard Python libraries
from os import path, listdir
import numpy as np
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(description="""
        RaMeDiES CH IND: provided a set of VCFs containing paternally and maternally inherited
        variants in probands, assess the per-individual significance of the deleteriousness of 
        the least expected compound heterozygous variant in a particular individual. The resulting 
        P-values in this case have to be corrected by the number of individuals analyzed rather 
        than by the number of genes. We recommend using this algorithm for smaller cohorts.
        For additional information, consult our paper and the manual.""")

    parser.add_argument('--variant_annots', type=str, default="CI", help="""
        String of codes for variant annotations:
        'C' for coding,
        'I' for intronic.
        Default: CI""")

    parser.add_argument('--i', type=str, help="""
        Input directory containing processed variant files. Should contain the slash symbol (/)
        at the end.""", default='')

    parser.add_argument('--o', type=str, default="out", help="""
        ID of the output files. Default: out""")

    parser.add_argument('--CADD_thr', type=float, default=1.5, help="""
        Minimal CADD score (raw score) for coding regions. Default: 1.5""")

    parser.add_argument('--SAI_thr', type=float, default=0.15, help="""
        Minimal SpliceAI score for coding regions. Default: 0.15""")

    parser.add_argument('--MAF', type=float, default=-1, help="""
        MAF threshold for variant filtering. Optional parameter for which
        the 'MAF' column in variant file has to be specified. Set to -1 for the absence of
        the MAF filter. default -1.""")

    # Boolean arguments
    parser.add_argument('--suppress_indels', help="""
        Include this argument to not count indel variants""", action='store_true')
    parser.set_defaults(suppress_indels=False)

    parser.add_argument('--no_qual_track', action='store_true', help="""
        Include if the input variant files are already quality-controlled and/or do not include
        the quality control column.""")
    parser.set_defaults(no_qual_track=False)

    return parser.parse_args()


def check_parameters(args):

    print('Checking command-line parameters...')
    if args.suppress_indels:
        print("> Indel suppression enabled (with --suppress_indels)")

    # Threshold output
    if 'C' in args.variant_annots:
        print(f"> Using raw CADD threshold of {args.CADD_thr}, change with --CADD_thr=<>")

    if 'I' in args.variant_annots:
        print(f"> Using SpliceAI threshold of {args.SAI_thr}, change with --SAI_thr=<>")
    print("")


if __name__ == "__main__":

    # Step 1: process and check parameters:
    args = parse_arguments()
    if not args.i.endswith('/'):
        args.i += '/'
    check_parameters(args)

    # Scores are stored in score_threshold_dict object throughout
    score_thr_dict = {'C': args.CADD_thr, 'I': args.SAI_thr}

    # List of variant annotations (SNP/indel, coding/intronic)
    # consequence_list and variant_annots are used interchangeably
    consequence_list = list(args.variant_annots)

    # Creating gene_instances object
    Gene_inst_dict = iol.parse_variant_scores_files(score_thr_dict,
                                                    consequence_list,
                                                    args.suppress_indels)

    # Total number of resulting genes. Should be around 16400
    N_genes = sl.count_genes(Gene_inst_dict)

    # Loading a list of pseudogenes, RNA genes and overlapping genes
    pseudogene_dict = iol.make_pseudogene_dict()

    # Total number of resulting genes. Should be around 16400
    ENS_ID_dict = sl.count_genes(Gene_inst_dict, return_dict=True)

    ENS2GeneID_dict = iol.make_ENS2GeneID_dict()

    # Initializing output file
    outfile_name = f"{args.o}_{cfg.CH_IND_result}.txt"
    with open(outfile_name, 'w') as outh:
        N_cohort = 0
        print(f"Reading input variant files from {args.i}")
        # Printing the output header
        sl.print_CH_IND_line(outh, first_line=True)
        for filename in listdir(args.i):
            # Processing each individual variant file
            # varcounts: inheritance ID (specified in cfg.inherited_from_dict) ->
            #	variant annotation (specified in cfg.var_annot_list) ->
            #	Number of per-individual variants
            varcounts = pvl.parse_variant_input(args.i + filename,
                                                Gene_inst_dict,
                                                maf_threshold=args.MAF,
                                                score_threshold_dict=score_thr_dict,
                                                pseudogene_dict=pseudogene_dict,
                                                consequence_list=consequence_list,
                                                suppress_indels_flag=args.suppress_indels,
                                                de_novo_flag=False,
                                                no_qual_track_flag=args.no_qual_track)

            # If no variants could be obtained from a input variant file, move to the next input file
            if varcounts != None:
                N_cohort += 1
            else:
                continue

            # muttarg_dict: mutational target -> CH_variant object (specified in stat_lib)
            # used for the inference of the least-expected CH per individual
            muttarg_dict = {}
            for ENS_ID in ENS_ID_dict.keys():
                # Initialization of VariantCollection objects
                # CH_variant objects are initialized by this command as well
                #	and stored in CH_list attribute of VariantCollection
                VariantCollection_obj = sl.VariantCollection(ENS_ID, Gene_inst_dict)
                VariantCollection_obj.make_CH_var_list(Gene_inst_dict)
                if VariantCollection_obj.CH_list == []:
                    continue

                muttarg = VariantCollection_obj.CH_list[0].mu
                muttarg_dict[muttarg] = VariantCollection_obj.CH_list[0]

            # Continue if no Ch variants were found
            if muttarg_dict == {}:
                # Delete all variant information from gene_instances
                Gene_inst_dict = sl.clean_Gene_inst_dict(Gene_inst_dict)
                continue

            top_mu = sorted(list(muttarg_dict.keys()))[0]
            # top_CH: CH_variant class (stat_lib) object
            #	corresponding to the least expected CH
            top_CH = muttarg_dict[top_mu]

            # Calculate statistic, write results to the output
            sl.calc_CH_IND_stat(outh,
                                varcounts,
                                top_CH,
                                Gene_inst_dict,
                                ENS2GeneID_dict)

            # Delete all variant information from gene_instances
            Gene_inst_dict = sl.clean_Gene_inst_dict(Gene_inst_dict)

    print(f"Read {N_cohort} processed variant files from {args.i}")
    print(f"Bonferroni P-value correction factor: {N_cohort}")
    print(f"Results written to: {outfile_name}")

    # Written on 01.23.2024 by Mikhail Moldovan, HMS DBMI
