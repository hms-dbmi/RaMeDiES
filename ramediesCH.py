"""
ramediesCH

Main script for the RaMeDiES cohort-level compound heterozygote analysis

Python version: 3.6 and above

Module requirements:
1. numpy
2. argparse
3. os

Next goes the description of the program. For more information about the usage,
consult the manual on GitHub and the help menu via the command:
python ramediesCH.py -h

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
Second, the VariantCollection objects (specified in stat_lib) as initiated 
for each gene. Unlike Gene objects, VariantCollection objects are not tied to
specific variant annotations, which allows to infer all CH variants by
using the make_comphet_var_list method of VariantCollection class and to calculate
the test statistic (called y) by the calc_comphet_y method of VariantCollection.
The products of per-cohort variant counts are also calculated.

Alternatively, under --metadata_run_mode enabled, the y values along with
the variant count products are loaded from the lists of intermediate outputs
the identifiers of which have to be provided by the user.

Parameters such as the gene count (num_genes) and the cohort size (N_probands)
are also calculated. N_probands should be specified in the matadata run 
through the --N_probands argument.

######
Third, if --metadata_write_mode is not enabled, will calculate statistics as 
described in the paper.

If you have any issues, contact us at mikhail_moldovan@hms.harvard.edu
"""

import init_objs_lib
import cfg
import parse_variant_lib
import stat_lib
from os import path, listdir
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="""RaMeDiES compound hets: provided a set of tab-delimited variant files containing 
        paternally and maternally inherited variants in probands, assess the significance of the enrichment 
        and deleteriousness of compound heterozygous variants landing within each gene. For additional information, 
        consult our paper and the GitHub Wiki page.""")

    parser.add_argument('--variant_annots', type=str, default="CI",
                        help="""String of codes for variant annotations:
        'C' for coding,
        'I' for intronic.
        Default: CI""")

    parser.add_argument('--i', type=str,
                        help="""Input directory containing processed variant files.""", default='')

    parser.add_argument('--M', type=str,
                        help="""Comma-separated list of prefixes of metadata files.
        needed only if metadata_run_mode is enabled""", default='')

    parser.add_argument('--o', type=str, default="out",
                        help="""Prefix of the output files. Default: out""")

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
    
    parser.add_argument('--SAI_thr', type=float, default=0.15,
        help="""Minimal SpliceAI score for splicing-affecting regions. 
        Default: 0.15""")

    parser.add_argument('--MAF', type=float, default=-1,
        help="""Maximal MAF threshold for variant filtering. Optional parameter for which
        the 'MAF' column in input variant files has to be specified. Set to -1 for the absence of
        the MAF filter. Default -1.""")
    
    parser.add_argument(
        '--coding_score', type=str, default="CADD",
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

    parser.add_argument('--N_probands', type=int, default=-1,
        help="""Specify the cohort size in the case of metadata_run_mode.
        The parameter is needed only for the false diagnosis rate estimates.
        If set to -1, false diagnosis rate will not be estimated. Default -1.""")

    # boolean arguments:
    parser.add_argument(
        '--missense_run', default=False, action='store_true',
        help="""Included if a missense-only predictor (e.g., PrimateAI-3D, AlphaMissense, REVEL) is
                specified via --coding_score=<>. Default: False""")

    parser.add_argument('--suppress_indels', default=False, action='store_true',
                        help="""Include this argument to not count indel variants. Default: False.""")

    parser.add_argument('--no_qual_track', default=False, action='store_true',
                        help="""Include if the input variant files are already quality-controlled 
        and/or do not include the quality control column. Default: False.""")

    parser.add_argument('--metadata_write_mode', default=False, action='store_true',
                        help="""Include this argument if only metadata files are needed. 
        In a default run, the program will produce metadata files anyway.
        Will produce an error if paired with --metadata_run_mode. Default: False.""")

    parser.add_argument('--metadata_run_mode', default=False, action='store_true',
                        help="""Include this argument for a metadata-only run.
        In a default run, a directory with variant files is specified in the --i parameter.
        In a metadata-only run, --i parameter is not specified and the input
        is given by the --M parameter. Default: False.""")

    parser.add_argument('--force_overwrite', default=False, action='store_true',
                        help="""Include to overwrite metadata files if the metadata files with 
        ID given by --o are present. If paired with --metadata_run_mode,
        this parameter will not affect anything. Default: False.""")

    return parser.parse_args()


def check_parameters(args):
    """
    :param args: object parsed using parse_args
    :return: NEW args object with modifications (if any were necessary)
    """

    print('Checking command-line parameters...')
    if args.suppress_indels:
        print("> Indel suppression enabled (with --suppress_indels)")

    if args.no_qual_track:
        print("> Quality track check disabled (with --no_qual_track)")

    if args.metadata_write_mode:
        print("> Metadata write mode enabled (with --metadata_write_mode)")
        print(f"  | Output will be printed to {args.o}_comphet_*.txt")

    if args.metadata_run_mode:
        print("> Metadata run enabled (with --metadata_run_mode)")
        print(f"  | Using the following files as input: {args.M}_*.txt")
        print(f"""  | Please make sure that the scores and thresholds used to produce the input 
                    match specifications of the run.""")

    if args.force_overwrite:
        print("> Force metadata overwrite enabled (with --force_overwrite)")

    if args.metadata_run_mode and args.M == '':
        raise AssertionError(
            "! Specify identifiers of metadata files using --M parameter (since --metadata_run_mode is enabled)")

    # metadata_write_mode and metadata_run_mode enabled, this gives an null run
    if args.metadata_write_mode and args.metadata_run_mode:
        raise AssertionError("! WARNING: NULL RUN because --metadata_write_mode and --metadata_run_mode are enabled.")

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

    # Checking the existence of intermediate metadata files.
    # If --force_overwrite is specified, those will be overwritten.
    varcount_file_exists = path.isfile(f"{args.o}_{cfg.mutnum_prod_CH}.txt")
    muttargs_file_exists = path.isfile(f"{args.o}_{cfg.muttargs_CH_ID}.txt")

    if varcount_file_exists and muttargs_file_exists:
        print(f"! Detected intermediate output file {args.o}_{cfg.mutnum_prod_CH}.txt")
        print(f"! Detected intermediate output file {args.o}_{cfg.muttargs_CH_ID}.txt")
        if args.force_overwrite:
            print("> Will overwrite intermediate outputs because --force_overwrite is specified")
        else:
            print(
                f"! Running in the metadata RUN mode using existing {args.o}_comphet_*.txt files. Rerun with --force_overwrite to overwrite instead.")
            args.metadata_run_mode = True
            args.M = args.o

    print("")
    return args


if __name__ == "__main__":

    # Step 1: parse command-line arguments
    args = parse_arguments()
    if args.i and not args.i.endswith('/'):
        args.i += '/'

    # Step 2: print specifications to terminal
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
    ENS_ID_dict = stat_lib.count_genes(Gene_inst_dict, return_dict=True)
    num_genes = len(ENS_ID_dict.keys())
    print("number of genes used in analysis:", num_genes)

    # Loading a list of pseudogenes, RNA genes and overlapping genes
    pseudogene_dict = init_objs_lib.make_pseudogene_dict()

    # Loading info from input variant files
    if not args.metadata_run_mode:

        N_cohort = 0  # compute this if we're in default or metadata_write_mode

        varcount_dict = {}
        print(f"Reading processed variant input files from {args.i}")
        for filename in listdir(args.i):

            varcounts = parse_variant_lib.parse_variant_input(args.i + filename,
                                                              Gene_inst_dict,
                                                              maf_threshold=args.MAF,
                                                              score_threshold_dict=score_thr_dict,
                                                              pseudogene_dict=pseudogene_dict,
                                                              consequence_list=consequence_list,
                                                              suppress_indels_flag=args.suppress_indels,
                                                              de_novo_flag=False,
                                                              no_qual_track_flag=args.no_qual_track,
                                                              missense_run_flag=args.missense_run)

            # If variant input file has been parsed sucessfully, respective proband is counted
            if varcounts:
                N_cohort += 1
            else:
                continue

            varcount_dict[filename] = varcounts

        print(f"Read {N_cohort} input variant files from {args.i}")

        # Writing distribution of by-annotation by-parental inheritance variant counts
        stat_lib.write_varcount_dist(varcount_dict, args.o, args.i, args_obj=args)

        # Calculating products of variant counts
        mutnum_prod_dict = stat_lib.mutnum_prod(varcount_dict, args.o, args_obj=args, args.i)

        # gene_comphet_mutdict: ENSEMBL ID -> VariantCollection object
        CH_dict = {}
        for ensembl_gene_id in ENS_ID_dict.keys():
            # Initializing VariantCollection objects
            VariantCollection_obj = stat_lib.VariantCollection(ensembl_gene_id, Gene_inst_dict)
            # Inference of CH variants
            VariantCollection_obj.make_comphet_var_list(Gene_inst_dict)
            # Inference of the main statistic (y)
            VariantCollection_obj.calc_comphet_y()
            # gene_comphet_mutdict update
            CH_dict[ensembl_gene_id] = VariantCollection_obj

        # Output CH mutational targets
        stat_lib.write_comphet_muttargs(CH_dict, args.o, args.i)

        # CH_y_dict: ensembl_gene_id -> CH y value
        CH_y_dict = stat_lib.calc_CH_y_dict(CH_dict)

    # Loading info from metadata files
    # Two objects must be loaded: y stats as CH_y_dict and products of single-annotation mut cnts as mutnum_prod_dict
    else:
        CH_dict = {}
        CH_y_dict = stat_lib.read_comphet_muttargs(args.M, args.variant_annots, args.suppress_indels)
        mutnum_prod_dict = stat_lib.read_mutnum_prods(args.M, args.variant_annots, args.suppress_indels)
        N_cohort = args.N_probands  # cohort size must be specified by user if we're not processing new files

    # actually compute the recurrence statistic (versus just producing the intermediate files)
    if not args.metadata_write_mode:
        stat_lib.calc_comphet_stat(CH_y_dict,
                                   CH_dict,
                                   Gene_inst_dict,
                                   mutnum_prod_dict,
                                   args.o,
                                   N_cohort)  # should have been computed in args.metadata_run_mode section

        print(f"Bonferroni P-value correction factor: {num_genes}")


# Written on 01.23.2024 by Mikhail Moldovan, HMS DBMI
