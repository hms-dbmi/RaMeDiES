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
import parse_variant_lib
import stat_lib
from os import path, listdir
import argparse


def parse_arguments():
    """
    :return: args object from parsing command-line arguments
    """

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="""RaMeDiES de novo: provided a set of tab-delimited variant files containing de novo variants in probands, 
        assess the significance of the enrichment and deleteriousness of variants landing 
        in genes hit by de novo mutations. For additional information, consult our paper and
        the GitHub Wiki page.""")
    
    parser.add_argument('--variant_annots', type=str, default="CI",
        help="""String of codes for variant annotations:
        'C' for coding,
        'I' for intronic.
        Default: CI""")
    
    parser.add_argument('--i', type=str,
        help="""Input directory containing processed variant files.""", default='')
    
    parser.add_argument('--M', type=str,
        help="""Comma-separated list of prefixes of metadata files.
        Needed only if metadata_run_mode is enabled""", default='')
    
    parser.add_argument('--o', type=str, default="out",
        help="""Prefix of the output files. Default: out""")
    
    parser.add_argument(
        '--coding_snv_thr', type=float, default=0.5,
        help="""Minimal deleteriousness score value for coding regions. Default: 0.5
        Suggested values for comphet variants by deleteriousness score type: 
        CADD (non Phred-scaled): 0.5
        AlphaMissense: 0.1
        PrimateAI-3D: 0.3
        REVEL: 0.5""")

    parser.add_argument(
        '--coding_indel_thr', type=float, default=0.5,
        help="""Minimal deleteriousness score for indels in coding regions.
        Non Phred-scaled CADD scores are always used for coding indels; this threshold may be different
        from the --coding_thr=<> specified if a deletriousness score other than CADD is selected via --coding_score=<>. 
        Default: 0.5 (value for raw CADD scores)""")

    parser.add_argument(
        '--SAI_thr', type=float, default=0.15,
        help="""Minimal SpliceAI score for intronic regions. 
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
    
    # Boolean arguments
    parser.add_argument('--missense_run', default=False, action='store_true',
        help="""Include this argument if PrimateAI-3D, AlphaMissense, REVEL, MisFitD or MisFitS is
                specified as coding_score.""", action='store_true')

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

    # de novo specific arguments:
    parser.add_argument('--gene_score', type=str, default="Genebayes",
        help="""Gene constraint score used in the weighted FDR procedure. 
        Values:
        1. 's_het_R': Roulette-corrected per-gene mean s_het estimates.
        2. 's_het_R_low95': Roulette-corrected per-gene estimates of 95 percent
            confidence interval lower bound for s_het.
        3. 's_het': mean s_het estimates from Cassa et al., 2017.
        4. 's_het_low95': estimates of 95 percent confidence interval lower bound for 
            s_het from Cassa et al., 2017.
        5. 'Genebayes': mean of s_het estimates from Zeng et al., 2023.
        6. 'Genebayes_low95': estimates of 95 percent confidence interval lower bound for 
            s_het from Zeng et al., 2023.
        7. 'LOEUF':  LOEUF estimates from Karczewski et al., 2020.
        8. 'PhyloP_prim': mean per-gene Primate PhyloP score.
        9. 'PhyloP_mamm': mean per-gene Mammal PhyloP score.
        10. 'MisFit_sgene_mis': per-gene MisFit missense estimates from Zhao et al., 2023.
        11. 'MisFit_sgene_ptv': per-gene MisFit PTV estimates from Zhao et al., 2023. 
        12. 'Mouse_knockout': gene stratification based on phenotypes of Mouse knockouts.
        13. 'OMIM_dom': only OMIM dominant genes are considered in the analyses. The effective
            number of tests is reduced by order of magnitude.
        Default: Genebayes""")

    parser.add_argument('--FDR', type=str, default="0.05,0.1",
        help="""Comma-separated FDR thresholds to flag and report findings. Default: 0.05,0.1""")

    parser.add_argument('--write_muttarg_sums', default=False, action='store_true',
        help="""Include to output per-gene sums of mutational targets instead of lists of
        per-variant mutational targets. Default: False.""")
    
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
        print(f"""  | Please make sure that the scores and thresholds used to produce the input 
                    match specifications of the run.""")

    if args.force_overwrite:
        print("> Force metadata overwrite enabled (with --force_overwrite)")

    if args.write_muttarg_sums:
        print("> Only sums of mutational targets will be reported (with --write_muttarg_sums)")

    if args.metadata_run_mode and args.M == '':
        raise AssertionError("! Specify identifiers of metadata files using --M parameter (since --metadata_run_mode is enabled)")

    # Input directory not specified or specified incorrectly
    # in the case of default run
    if not args.metadata_run_mode and (args.i == '' or not args.i.endswith('/')):
        raise AssertionError("! Input directory specified incorrectly: --i="+args.i)

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

    # Gene constraint metric output
    print(f"> Using {args.gene_score} as gene constraint metric in weighted FDR.")

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
            print(f"""! Running in the metadata RUN mode using existing {args.o}_denovo_*.txt files. 
                        Rerun with --force_overwrite to overwrite instead.""")
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
    score_thr_dict = {'C' : args.coding_snv_thr,
                      'I' : args.SAI_thr, 
                      'CInd' : args.coding_indel_thr}

    # List of variant annotations (SNP/indel, coding/intronic)
    # consequence_list and variant_annots are used interchangeably
    consequence_list = list(args.variant_annots)

    # Creating gene_instances object
    gene_instances_dict = init_objs_lib.parse_variant_scores_files(score_thr_dict,
                                                                   consequence_list,
                                                                   args.coding_score,
                                                                   args.suppress_indels)

    # Loading a list of pseudogenes, RNA genes and overlapping genes
    pseudogene_dict = init_objs_lib.make_pseudogene_dict()

    # Total number of resulting genes. Should be around 16400
    num_genes = stat_lib.count_genes(gene_instances_dict)
    print("number of genes used in analysis:", num_genes)

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
                                                              no_qual_track_flag=args.no_qual_track,
                                                              missense_run_flag=args.missense_run)

            # If variant file has been parsed successfully, respective proband is counted
            if varcounts is not None:
                num_samples += 1
            else:
                continue

            varcount_dict[filename] = varcounts

        print(f"Read {num_samples} processed variant input files from {args.i}")

        # Initializing by_annot_varcount_dict
        # Writing metadata and QC files
        by_annot_varcount_dict = stat_lib.make_by_annot_varcount_dict(varcount_dict, args.o, args_obj = args)
        stat_lib.write_varcount_dist(varcount_dict, args.o, args.i, args_obj = args)
        stat_lib.write_muttargs(gene_instances_dict, args.o, args.i, args_obj = args)

    # Loading info from metadata files
    else:
        gene_instances_dict = stat_lib.load_muttargs_from_filelist(gene_instances_dict, args.M, args.variant_annots,
                                                                   args.suppress_indels)
        by_annot_varcount_dict = stat_lib.load_by_annot_varcount_dict(args.M, args.variant_annots, args.suppress_indels)

        # N_probands should be specified by the user in the case of metadata_run_mode
        num_samples = args.N_probands

    # Calculating statistics
    if not args.metadata_write_mode:
        # ensg_to_genename: ensembl_gene_id -> list of Variant objects
        gene_mutdict = stat_lib.make_gene_mutdict(gene_instances_dict)

        # ENS_ID_y_dict: ensembl_gene_id -> y statistic
        ENS_ID_y_dict = stat_lib.count_y(gene_instances_dict)

        # Main function for the statistical inference. Look in stat_lib for details
        stat_lib.calc_denovo_stat(ENS_ID_y_dict,
                                  by_annot_varcount_dict,
                                  gene_instances_dict,
                                  args.o,
                                  num_samples,
                                  gene_mutdict,
                                  args.FDR,
                                  args.gene_score,
                                  num_genes)

        print(f"Bonferroni P-value correction factor: {num_genes}")

# Written on 01.23.2024 by Mikhail Moldovan, HMS DBMI
