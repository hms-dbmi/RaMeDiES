"""
ramediesDN

Prerequisites:
* Python version: 3.6 and above
* Python packages: numpy, argparse, os

Steps:
1. Create a collection of "gene objects" (see stat_lib) stratified by variant type (coding/intronic, SNV/indel) and
   by Ensembl gene ID
2. Populate these gene objects with observed de novo variants in one of two ways:
   (a) read in variants from the specified processed input files
   (b) read in variants from a "by_annot_varcount_dict" object (created by "make_by_annot_varcount_dict") in stat_lib
3. Print gene count (num_genes), cohort size (N_probands) to terminal
4. Create output files as specified in the user manual:
   (a) Default and metadata_write_mode will export a variant_counts, variant_distribution, and mutational_targets file
   (b) Default and metadata_run_mode will export a cohort_recurrence file

If you have any issues, contact us at mikhail_moldovan@hms.harvard.edu
"""

import init_objs_lib
import cfg
import parse_VCF_lib as parse_variant_lib
import stat_lib
from os import path, listdir
import argparse


def parse_arguments():
    """
    :return: args object from parsing command-line arguments
    """

    parser = argparse.ArgumentParser(description="""
        RaMeDiES de novo: provided a set of VCFs containing de novo variants in probands, 
        assess the significance of the enrichment and deleteriousness of variants landing 
        in genes hit by de novo mutations. For additional information, consult our paper and
        the manual.""")

    parser.add_argument('--variant_annots', type=str, default="CI", help="""
        String of codes for variant annotations:
        'C' for coding,
        'I' for intronic.
        Default: CI""")

    parser.add_argument('--i', type=str, help="""
        Input directory containing processed variant files.""", default='')

    parser.add_argument('--M', type=str, help="""
        Comma-separated list of prefixes of metadata files.
        needed only if metadata_run_mode is enabled""", default='')

    parser.add_argument('--o', type=str, default="out", help="""
        ID of the output files. Default: out""")

    parser.add_argument('--CADD_thr', type=float, default=0.5, help="""
        Minimal CADD score (raw score) for coding regions. Default: 0.5""")

    parser.add_argument('--SAI_thr', type=float, default=0.05, help="""
        Minimal SpliceAI score for coding regions. Default: 0.05""")

    parser.add_argument('--FDR', type=str, default="0.05,0.1", help="""
        Comma-separated FDR values. Default: 0.05,0.1""")

    parser.add_argument('--MAF', type=float, default=-1, help="""
        MAF threshold for variant filtering. Optional parameter for which
        the 'MAF' column in input variant files has to be specified. Set to -1 for the absence of
        the MAF filter. default -1.""")

    parser.add_argument('--N_probands', type=int, default=-1, help="""
        Specify the cohort size in the case of metadata_run_mode.
        The parameter is needed only for the false diagnosis rate estimates.
        If set to -1, false diagnosis rate will just not be estimated. Default -1.""")

    # Boolean arguments
    parser.add_argument('--suppress_indels', help="""
        Include this argument to not count indel variants""", action='store_true')
    parser.set_defaults(suppress_indels=False)

    parser.add_argument('--metadata_write_mode', action='store_true', help="""
        Include this argument if only metadata files are needed. 
        In a default run, the program will produce metadata files anyway.
        Will produce an error if paired with --metadata_run_mode.""")
    parser.set_defaults(metadata_write_mode=False)

    parser.add_argument('--metadata_run_mode', action='store_true', help="""
        Include this argument for a metadata-only run.
        In a default run, a directory with variant files is specified in the --i parameter.
        In a metadata-only run, --i parameter is not specified and the input
        is given by the --M parameter.""")
    parser.set_defaults(metadata_run_mode=False)

    parser.add_argument('--force_overwrite', action='store_true', help="""
        Include to overwrite metadata files if the metadata files with 
        ID given by --o are present. If paired with --metadata_run_mode,
        this parameter will not affect anything.""")
    parser.set_defaults(force_overwrite=False)

    parser.add_argument('--no_qual_track', action='store_true', help="""
        Include if the input variant files are already quality-controlled and/or do not include
        the quality control column.""")
    parser.set_defaults(no_qual_track=False)

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
        print(f"  | Output will be printed to {args.o}_denovo_*.txt")

    if args.metadata_run_mode:
        print("> Metadata run enabled (with --metadata_run_mode)")
        print(f"  | Using the following files as input: {args.M}_*.txt")

    if args.force_overwrite:
        print("> Force metadata overwrite enabled (with --force_overwrite)")

    if args.metadata_run_mode and args.M == '':
        raise AssertionError("! Specify identifiers of metadata files using --M parameter (since --metadata_run_mode is enabled)")

    # Input directory not specified or specified incorrectly
    # in the case of default run
    if not args.metadata_run_mode and (args.i == '' or not args.i.endswith('/')):
        raise AssertionError("! Input directory specified incorrectly: --i="+args.i)

    # metadata_write_mode and metadata_run_mode enabled, this gives an null run
    if args.metadata_write_mode and args.metadata_run_mode:
        raise AssertionError("! WARNING: NULL RUN because --metadata_write_mode and --metadata_run_mode are enabled.")

    # Threshold output
    if 'C' in args.variant_annots:
        print(f"> Using raw CADD threshold of {args.CADD_thr}, change with --CADD_thr=<>")

    if 'I' in args.variant_annots:
        print(f"> Using SpliceAI threshold of {args.SAI_thr}, change with --SAI_thr=<>")

    # Checking the existence of intermediate metadata files.
    # If --force_overwrite is specified, those will be overwritten.
    varcount_file_exists = path.isfile(f"{args.o}_{cfg.varcount_sums_DN}.txt")
    muttargs_file_exists = path.isfile(f"{args.o}_{cfg.muttargs_list_DN_ID}.txt")

    if varcount_file_exists and muttargs_file_exists:
        print(f"! Detected intermediate output file {args.o}_{cfg.varcount_sums_DN}.txt")
        print(f"! Detected intermediate output file {args.o}_{cfg.muttargs_list_DN_ID}.txt")
        if args.force_overwrite:
            print("> Will overwrite intermediate outputs because --force_overwrite is specified")
        else:
            print(f"! Running in the metadata RUN mode using existing {args.o}_denovo_*.txt files. Rerun with --force_overwrite to overwrite instead.")
            args.metadata_run_mode = True
            args.M = args.o
    print("")

    return args


if __name__ == "__main__":

    # Step 1: parse command-line arguments
    args = parse_arguments()
    if not args.i.endswith('/'):
        args.i += '/'

    # Step 2: print specifications to terminal
    args = check_parameters(args)

    # Scores are stored in score_threshold_dict object throughout
    score_thr_dict = {'C': args.CADD_thr, 'I': args.SAI_thr}

    # List of variant annotations (SNP/indel, coding/intronic)
    # consequence_list and variant_annots are used interchangeably
    consequence_list = list(args.variant_annots)

    # Creating gene_instances object
    gene_instances_dict = init_objs_lib.parse_variant_scores_files(score_thr_dict,
                                                                   consequence_list,
                                                                   args.suppress_indels)

    # Loading a list of pseudogenes, RNA genes and overlapping genes
    pseudogene_dict = init_objs_lib.make_pseudogene_dict()

    # Total number of resulting genes. Should be around 16400
    num_genes = stat_lib.count_genes(gene_instances_dict)
    print("number of genes:", num_genes)

    if not args.metadata_run_mode:
        varcount_dict = {}
        num_samples = 0
        print(f"Reading processed variant input files from {args.i}")
        for filename in listdir(args.i):

            varcounts = parse_variant_lib.parse_variant_input(args.i + filename,
                                                              gene_instances_dict,
                                                              maf_threshold=args.MAF,
                                                              score_threshold_dict=score_thr_dict,
                                                              pseudogene_dict=pseudogene_dict,
                                                              consequence_list=consequence_list,
                                                              suppress_indels_flag=args.suppress_indels,
                                                              de_novo_flag=True,
                                                              no_qual_track_flag=args.no_qual_track)

            # If variant file has been parsed successfully, respective proband is counted
            if varcounts is not None:
                num_samples += 1
            else:
                continue

            varcount_dict[filename] = varcounts

        print(f"Read {num_samples} processed variant input files from {args.i}")

        # Initializing by_annot_varcount_dict
        # Writing metadata and QC files
        by_annot_varcount_dict = stat_lib.make_by_annot_varcount_dict(varcount_dict, args.o)
        stat_lib.write_varcount_dist(varcount_dict, args.o, args.variant_annots, score_thr_dict, args.i)
        stat_lib.write_muttargs(gene_instances_dict, args.o, args.i)

    # Loading info from metadata files
    else:
        gene_instances_dict = stat_lib.load_muttargs_from_filelist(gene_instances_dict, args.M, args.variant_annots,
                                                                   args.suppress_indels)
        by_annot_varcount_dict = stat_lib.load_by_annot_varcount_dict(args.M, args.variant_annots, args.suppress_indels)

        # N_probands should be specified by the user in the case of metadata_run_mode
        num_samples = args.N_probands

    # Calculating statistics
    if not args.metadata_write_mode:
        # gene_mutdict: ensembl_gene_id -> list of Variant objects
        gene_mutdict = stat_lib.make_gene_mutdict(gene_instances_dict)

        # ENS_ID_y_dict: ensembl_gene_id -> y statistic
        ENS_ID_y_dict = stat_lib.count_y(gene_instances_dict, args.o)

        # Main function for the statistical inference. Look in stat_lib for details
        stat_lib.calc_DN_stat(ENS_ID_y_dict,
                              by_annot_varcount_dict,
                              gene_instances_dict,
                              args.o,
                              num_samples,
                              gene_mutdict,
                              args.FDR,
                              num_genes)

        print(f"Bonferroni P-value correction factor: {num_genes}")

# Written on 01.23.2024 by Mikhail Moldovan, HMS DBMI
