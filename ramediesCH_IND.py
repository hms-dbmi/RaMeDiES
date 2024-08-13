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
using the make_comphet_var_list method of VariantCollection class and to calculate
the test statistic (called y) by the calc_comphet_y method of VariantCollection.
The products of per-cohort variant counts are also calculated. Finally,
a statistic described in our paper is calculated and the result is written to the
output file with identified provided through the --o option.

If you have any issues, contact us at mikhail_moldovan@hms.harvard.edu
"""

## Libraries 
# import init_objs_lib as iol
# import cfg
# import parse_VCF_lib as pvl
# import stat_lib as sl
#
## Standard Python libraries
# from os import path, listdir
# import numpy as np
# import argparse


import init_objs_lib
import cfg
import parse_variant_lib
import stat_lib
from os import path, listdir
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="""RaMeDiES CH IND: provided a set of tab-delimited variant files containing paternally 
        and maternally inherited variants in probands, assess the per-individual significance of the 
        deleteriousness of the single least expected compound heterozygous variant in a particular individual. 
        The resulting P-values in this case have to be corrected by the number of individuals analyzed rather 
        than by the number of genes. We recommend using this algorithm for smaller cohorts. For additional 
        information, consult our paper and the GitHub Wiki page.""")

    parser.add_argument(
        '--variant_annots', type=str, default="CI", choices=['C', 'I', 'CI', 'IC'], help="""
        String of codes for variant annotations:
        'C' for coding,
        'I' for intronic.
        Default: CI""")

    parser.add_argument(
        '--i', type=str,
        help="""Input directory containing processed input variant files.""", default='')

    parser.add_argument(
        '--o', type=str, default="CHIND_result",
        help="""Prefix for output files""")

    parser.add_argument(
        '--coding_snv_thr', type=float, default=1.5,
        help="""Minimal deleteriousness score value for coding regions. Default: 1.5
        Suggested values for comphet variants by deleteriousness score type: 
        CADD (non Phred-scaled): 1.5
        AlphaMissense: 0.1
        PrimateAI-3D: 0.3
        REVEL: 0.5""")

    parser.add_argument(
        '--coding_indel_thr', type=float, default=1.5,
        help="""Minimal deleteriousness score for indels in coding regions.
        Non Phred-scaled CADD scores are always used for coding indels; this threshold may be different
        from the --coding_thr=<> specified if a deletriousness score other than CADD is selected via --coding_score=<>. 
        Default: 1.5 (value for raw CADD scores)""")

    parser.add_argument(
        '--SAI_thr', type=float, default=0.15,
        help="""Minimal SpliceAI score for intronic regions. 
        Default: 0.15""")

    parser.add_argument(
        '--MAF', type=float, default=-1,
        help="""Maximal MAF threshold for variant filtering. Optional parameter for which
        the 'MAF' column in input variant files has to be specified. Set to -1 for the absence of
        the MAF filter. Default -1.""")

    parser.add_argument(
        '--coding_score', type=str, default="CADD", choices=['CADD', 'AlphaMissense', 'PAI3D', 'REVEL'],
        help="""Deleteriousness score type for coding SNP variants.
        Note: non-Phred-scaled CADD scores are ALWAYS used for scoring coding indels 
              (unless the --suppress_indels flag is indicated, in which case no coding indels are included.)
        Note: Some deleteriousness scores (marked with an *) only score missense SNVs (see --missense_run).
        Possible values for --coding_score=<>:
        - CADD : raw CADD score
        - AlphaMissense *: AlphaMissense score
        - PAI3D *: Variant deleteriousness predictions from PrimateAI-3D.
        - REVEL *: REVEL score
        Default: CADD""")

    # Boolean arguments
    parser.add_argument(
        '--missense_run', default=False, action='store_true',
        help="""Included if a missense-only predictor (e.g., PrimateAI-3D, AlphaMissense, REVEL) is
                specified via --coding_score=<>. Default: False""")

    parser.add_argument('--suppress_indels', default=False, action='store_true',
                        help="""Include this argument to not count indel variants. Default: False.""")

    parser.add_argument('--no_qual_track', default=False, action='store_true',
                        help="""Include if the input variant files are already quality-controlled 
        and/or do not include the quality control column. Default: False.""")

    return parser.parse_args()


def check_parameters(args):
    """
    :param args: dictionary of command-line arguments parsed by parser
    :return:
    """

    print('Checking command-line parameters...')
    if args.suppress_indels:
        print("> Indel suppression enabled (with --suppress_indels)")

    if args.no_qual_track:
        print("> Quality track check disabled (with --no_qual_track)")

    # Score restrictions need to match the rules of processing of variant collections.
    if 'C' in args.variant_annots:
        print(f"> Using {args.coding_score} to process coding SNVs with threshold {args.coding_snv_thr}.")
        print(f"| Change score type with --coding_score=<> and threshold with --coding_snv_thr=<>.")

        if args.missense_run:
            print("| --missense_run enabled, only SNVs with a missense impact are considered")

        if args.coding_score != "CADD":
            print(f"| Consider updating the \"coding_snv_score\" key/value pair in \"vcf_format_dict\" in cfg ")
            print(f"  to read in the correct columm (i.e., not \"CADD-raw\") from your input variant files.")

            if args.coding_score in ['PAI3D', 'AlphaMissense', 'REVEL'] and not args.missense_run:
                print(
                    f"! WARNING: {args.coding_score} only scores missense SNVs. NO OTHER coding SNVs will be considered.")
                args.missense_run = True

            # No deleteriousness scores besides CADD are available for short indels
            if not args.suppress_indels:
                print(f"| Using CADD for coding indels with non Phred-scaled threshold {args.coding_indel_thr}.")
                print(" ! Change coding indel score threshold with --coding_indel_thr=<>.")
                print(" ! Include --suppress_indels flag to consider only SNVs in coding and intronic regions.")

        elif not args.suppress_indels and args.coding_snv_thr != args.coding_indel_thr:
            print(
                f"! NOTE: The CADD threshold for SNVs ({args.coding_snv_thr}) is not the same CADD threshold for indels ({args.coding_indel_thr})")
            print("  Change these value(s) with --coding_snv_thr=<> or --coding_indel_thr=<>")

    if 'I' in args.variant_annots:
        print(f"> Using SpliceAI to process intronic SNVs/indels with threshold {args.SAI_thr}.")
        print(f"| Change intronic variant score threshold with --SAI_thr=<>.")
        print(f"| Exclude all intronic variants with --variant_annots=C.")

    print(f"")
    return args


if __name__ == "__main__":

    # Step 1: process and check parameters:
    args = parse_arguments()
    if not args.i.endswith('/'):
        args.i += '/'
    args = check_parameters(args)

    # Scores are stored in score_threshold_dict object throughout
    score_thr_dict = {'C': args.coding_snv_thr,
                      'I': args.SAI_thr,
                      'CInd': args.coding_indel_thr}

    # List of variant annotations (SNP/indel, coding/intronic)
    # consequence_list and variant_annots are used interchangeably
    consequence_list = list(args.variant_annots)

    # Creating gene_instances object
    Gene_inst_dict = init_objs_lib.parse_variant_scores_files(score_thr_dict,
                                                              consequence_list,
                                                              args.coding_score,
                                                              args.suppress_indels)

    # Total number of resulting genes. Should be around 16400
    num_genes = stat_lib.count_genes(Gene_inst_dict)
    print("number of genes used in analysis:", num_genes)

    # Loading a list of pseudogenes, RNA genes and overlapping genes
    pseudogene_dict = init_objs_lib.make_pseudogene_dict()

    # Total number of resulting genes. Should be around 16400
    genes_to_include = stat_lib.count_genes(Gene_inst_dict, return_dict=True)

    ensg_to_genename = init_objs_lib.make_ENS2GeneID_dict()

    # Initializing output file
    outfile_name = f"{args.o}_{cfg.CH_IND_result}.txt"
    with open(outfile_name, 'w') as out_handle:

        N_cohort = 0
        print(f"Reading input variant files from {args.i}")

        # Printing the output header
        stat_lib.print_comphet_individual_line(out_handle, first_line=True)

        # Processing each individual variant file
        for filename in listdir(args.i):

            # variant_counts: inheritance ID (specified in cfg.inherited_from_dict) ->
            #	variant annotation type (specified in cfg.var_annot_list) ->
            #	Number of per-individual variants
            variant_counts = parse_variant_lib.parse_variant_input(args.i + filename,
                                                                   Gene_inst_dict,
                                                                   maf_threshold=args.MAF,
                                                                   score_threshold_dict=score_thr_dict,
                                                                   pseudogene_dict=pseudogene_dict,
                                                                   consequence_list=consequence_list,
                                                                   suppress_indels_flag=args.suppress_indels,
                                                                   de_novo_flag=False,
                                                                   no_qual_track_flag=args.no_qual_track,
                                                                   missense_run_flag=args.missense_run)

            # If no variants could be obtained from a input variant file, move to the next input file
            if variant_counts:
                N_cohort += 1
            else:
                continue

            # muttarg_dict: mutational target -> CH_variant object (specified in stat_lib)
            # used for the inference of the least-expected CH per individual
            muttarg_dict = {}

            # Initialization of VariantCollection objects
            for ensembl_gene_id in genes_to_include.keys():

                # CH_variant objects are initialized by this command as well
                #	and stored in CH_list attribute of VariantCollection
                VariantCollection_obj = stat_lib.VariantCollection(ensembl_gene_id, Gene_inst_dict)
                VariantCollection_obj.make_comphet_var_list(Gene_inst_dict)
                if VariantCollection_obj.CH_list == []:
                    continue

                muttarg = VariantCollection_obj.CH_list[0].mu
                muttarg_dict[muttarg] = VariantCollection_obj.CH_list[0]

            # Continue if no CH variants were found
            if muttarg_dict == {}:
                # Delete all variant information from gene_instances
                Gene_inst_dict = stat_lib.clean_Gene_inst_dict(Gene_inst_dict)
                continue

            top_mu = sorted(list(muttarg_dict.keys()))[0]
            # top_comphet: CH_variant class (stat_lib) object
            #	corresponding to the least expected CH
            top_comphet = muttarg_dict[top_mu]

            # Calculate statistic, write results to the output
            stat_lib.calc_comphet_individual_statistic(out_handle,
                                                       variant_counts,
                                                       top_comphet,
                                                       Gene_inst_dict,
                                                       ensg_to_genename)

            # Delete all variant information from gene_instances
            Gene_inst_dict = stat_lib.clean_Gene_inst_dict(Gene_inst_dict)

    print(f"Read {N_cohort} processed variant files from {args.i}")
    print(f"Bonferroni P-value correction factor: {N_cohort}")
    print(f"Results written to: {outfile_name}")

# Written on 01.23.2024 by Mikhail Moldovan, HMS DBMI
