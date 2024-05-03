"""
Functions used to run the pathway analysis presented in our paper.
Full instructions can be found in our wiki:
https://github.com/hms-dbmi/RaMeDiES/wiki/Pathway-analysis

Last updated: 2024-04-30
"""

import sys
import requests
import argparse
import cfg
import rpy2.robjects as robjects
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage
from phrank import Phrank


########################################################################################################

def create_phrank_instance(
    hpo_to_disease=cfg.script_directory+'/data/phrank/disease_to_pheno.build127.txt',
    gene_to_disease=cfg.script_directory+'/data/phrank/gene_to_disease.build127.txt',
    hpo_to_gene=cfg.script_directory+'/data/phrank/gene_to_pheno.amelie.txt',
    hpo_dag=cfg.script_directory+'/data/phrank/hpodag.txt'):
    """
    :param hpo_to_disease: tab-delimited HPO_ID \t disease_ID, downloaded  https://bitbucket.org/bejerano/phrank/src/
        master/demo/data/disease_to_pheno.build127.txt
    :param gene_to_disease: tab-delimited gene_id \t disease_id, downloaded from https://bitbucket.org/bejerano/
        phrank/src/master/demo/data/gene_to_disease.build127.txt
    :param hpo_to_gene: tab-delimited hpo_id \t gene_id, downloaded from https://bitbucket.org/bejerano/phrank/src/
        master/demo/data/gene_to_pheno.amelie.txt
    :param hpo_dag: HPO ontology represented as ordered pairs of child/parent, downloaded from
        https://bitbucket.org/bejerano/phrank/src/master/demo/data/hpodag.txt
    :return: a Phrank instance that can be used for multiple queries
    """

    return Phrank(hpo_dag,
                  diseaseannotationsfile=hpo_to_disease,
                  diseasegenefile=gene_to_disease,
                  geneannotationsfile=hpo_to_gene)


########################################################################################################

def patient_pheno_match(patient1_hpos, patient2_hpos, p_hpo):
    """
    :param patient1_hpos: set of HPO terms belonging to patient 1
    :param patient2_hpos: set of HPO terms belonging to patient 2
    :param p_hpo: Phrank instance
    :return: 3 scores-- the pairwise similarity, the max possible for patient 1, and the max for patient 2
    """

    patient1_max = p_hpo.compute_maximal_match(set([a.split('|')[0] for a in patient1_hpos]))
    patient2_max = p_hpo.compute_maximal_match(set([a.split('|')[0] for a in patient2_hpos]))
    pairwise_match = p_hpo.compute_phenotype_match(set([a.split('|')[0] for a in patient1_hpos]),
                                                   set([a.split('|')[0] for a in patient2_hpos]))

    return pairwise_match, patient1_max, patient2_max


########################################################################################################

def hpo_from_hpofile(hpo_file):
    """
    :param hpo_file: tab-delimited file of HPO terms extracted for UDN patients
    :return: dictionary of individual ID -> set of HPO term IDs (positive)
    """

    hpo_terms = {}

    with open(hpo_file) as hpo_handle:
        header = None
        for l in hpo_handle:
            if l.startswith('#'):
                continue
            elif not header:
                header = l.strip().split('\t')

                if 'patient_id' not in header:
                    sys.stderr.write('[ERROR] Could not find "patient_id" in HPO terms file: '+hpo_file+' ...exiting\n')
                    sys.exit(1)
                if 'hpo_id' not in header:
                    sys.stderr.write('[ERROR] Could not find "hpo_id" in HPO terms file: ' + hpo_file + ' ...exiting\n')
                    sys.exit(1)
                if 'hpo_name' not in header:
                    sys.stderr.write('Could not find "hpo_name" in HPO terms file: ' + hpo_file + ' ...listing term names as "unspecified"\n')
                    sys.exit(1)
                if 'patient_status' not in header:
                    sys.stderr.write('Could not find "patient_status" in HPO terms file: '+hpo_file+' ...assuming all patients are affected\n')
                if 'hpo_status' not in header:
                    sys.stderr.write('Could not find "hpo_status" in HPO terms file: ' + hpo_file + ' ...assuming all HPO terms are present\n')
                if 'quality_flag' not in header:
                    sys.stderr.write('Could not find "quality_flag" in HPO terms file: ' + hpo_file + ' ...assuming all HPO terms PASS\n')
                continue
            v = l[:-1].split('\t')

            if 'hpo_status' in header and v[header.index('hpo_status')].lower() not in ['true', 'yes', 'y', 'present', 'positive']:
                continue
            if 'patient_status' in header and v[header.index('patient_status')].lower() not in ['true', 'yes', 'y', 'affected']:
                continue
            if 'quality_flag' in header and v[header.index('quality_flag')].lower() not in ['true', 'yes', 'y', 'pass']:
                continue

            hpo_id = v[header.index('hpo_id')]
            hpo_name = 'unspecified'
            if 'hpo_name' in header:
                hpo_name = '-'.join(v[header.index('hpo_name')].replace(',', '').replace('|', '').split())

            sample_id = v[header.index('patient_id')]
            if sample_id not in hpo_terms:
                hpo_terms[sample_id] = set()
            hpo_terms[sample_id].add(hpo_id+'|'+hpo_name)

    return {individual_id: sorted(list(hpo_ids)) for individual_id, hpo_ids in hpo_terms.items()}


########################################################################################################

def pairwise_phrank_similarity(
    hpo_file=cfg.script_directory+'/test/phrank/test_patient_hpoterms.tsv',
    out_file=cfg.script_directory+'/test/phrank/test_pairwise_similarity.tsv'):
    """
    :param hpo_file: tab-delimited file with patient_id, affected_status, hpo_id, hpo_name, status (True or False), quality_flag (PASS or False)
    :param out_file: full path to a file where tab-delimited pairwise similarity scores should be written
    :return: (None) ALL-AGAINST-ALL phenotype comparisons of all patients with 1+ positive HPO terms from hpo_file
    """

    # how well does this gene match the patient's phenotypes?
    phrank_instance = create_phrank_instance(
        cfg.script_directory+'/data/phrank/pheno_to_disease_2021-11-22.txt',
        cfg.script_directory+'/data/phrank/gene_to_disease_2021-11-22.txt',
        cfg.script_directory+'/data/phrank/phenotype_to_gene_2021-11-22.txt',
        cfg.script_directory+'/data/phrank/hpo_dag_2021-11-22.txt')

    # first, get all the HPO terms per affected individual
    patient_to_hpos = hpo_from_hpofile(hpo_file)  # dictionary of patient_id -> {hpo_id1, hpo_id2, ...}

    sample_ids = sorted([udn_id for udn_id, hpo_set in patient_to_hpos.items() if len(hpo_set) > 0])

    outhandle = open(out_file, 'w')
    outhandle.write('\t'.join(['normalized_pairwise_similarity', 'raw_pairwise_similarity',
                               'sorted_patient_max', 'patient1', 'patient2'])+'\n')

    for i in range(len(sample_ids)):
        for j in range(i, len(sample_ids)):
            if i < j:
                pairwise, patient1max, patient2max = patient_pheno_match(
                    patient_to_hpos.get(sample_ids[i], set()),
                    patient_to_hpos.get(sample_ids[j], set()),
                    phrank_instance)
                outhandle.write('\t'.join([
                    str(pairwise/(max(patient1max, patient2max))),
                    str(pairwise),
                    str(max(patient1max, patient2max))+','+str(min(patient1max, patient2max)),
                    sample_ids[i],
                    sample_ids[j]
                ])+'\n')
    outhandle.close()


########################################################################################################

def create_patient_clusters(pairwise_similarities=cfg.script_directory+'/test/phrank/test_pairwise_similarity.tsv',
                            cluster_assignments=cfg.script_directory+'/test/phrank/test_cluster_assignments.tsv'):
    """
    :param pairwise_similarities: tab-delimited file created in "pairwise_phrank_similarity"
    :param cluster_assignments: tab-delimited output file with sample->cluster mapping
    :return: (None) assignments of samples to clusters based on complete-linkage hierarchical clustering and cutting the
        dendrogram at height 3.5
    """

    # Code written in R to perform hierarchical clustering
    r_code = """
    create_clusters <- function(infile, outfile, linkage_type="complete", cut_height=3.5) {
        library(cluster)
        
        # Read in Phrank similarity scores produced in previous step
        data <- read.table(infile, header = TRUE, sep = "\t")
        items <- unique(c(data$patient1,data$patient2))

        # Initialize an empty similarity matrix
        n_items <- length(items)
        similarity_matrix <- matrix(1, nrow = n_items, ncol = n_items)

        # Fill in the similarity matrix
        for (i in 1:nrow(data)) {
            item1 <- data[i, "patient1"]
            item2 <- data[i, "patient2"]
            similarity <- data[i, "normalized_pairwise_similarity"]
            
            # Find the row and column indices for the items
            row_index <- which(items == item1)
            col_index <- which(items == item2)
  
            # Fill in the similarity values in the matrix
            similarity_matrix[row_index, col_index] <- similarity
            similarity_matrix[col_index, row_index] <- similarity  # Similarity matrix is symmetric
        }
        rownames(similarity_matrix) <- colnames(similarity_matrix) <- items

        # Calculate the distance matrix from the similarity matrix
        distance_matrix <- 1 - similarity_matrix

        # Perform hierarchical clustering
        clustering_result <- agnes(distance_matrix, method = linkage_type)
        clustering_result <- as.hclust(clustering_result)  # Ensure it's in hclust format

        # Cut the dendrogram to obtain cluster assignments and write to file
        cluster_assignments <- cutree(clustering_result, h = cut_height)
        write.table(cluster_assignments, 
            file = outfile, 
            sep = "\t", quote = FALSE, col.names=FALSE)
    }
    """

    # Define then call and execute the clustering function
    robjects.r(r_code)
    robjects.r['create_clusters'](pairwise_similarities, cluster_assignments)


########################################################################################################

def create_gene_groups(candidate_genes_list=cfg.script_directory+'/test/phrank/test_candidate_genes.tsv',
                       cluster_assignments=cfg.script_directory+'/test/phrank/test_cluster_assignments.tsv',
                       cluster_details=cfg.script_directory+'/test/phrank/test_cluster_details.tsv'):
    """
    :param candidate_genes_list: tab-delimited file containing patient IDs and their candidate genes (to be used in
        gene set enrichment analysis) with other relevant information as desired
    :param cluster_assignments: tab-delimited file produced by "create_patient_clusters" function assigning patients
        to clusters based on the semantic similarity of their phenotype terms
    :param cluster_details: output tab-delimited file listing cluster_id, sample_id, gene_name to be used as a query
        for gene set enrichment analysis
    :return: (None) tab-delimited file with clusters -> genes to run gene set enrichment
    """

    # Read cluster assignments
    cluster_to_patients = {}
    with open(cluster_assignments) as chandle:
        for cline in chandle:
            sample_id, cluster_id = cline.strip().split()
            if 'HC'+cluster_id not in cluster_to_patients:
                cluster_to_patients['HC'+cluster_id] = set()
            cluster_to_patients['HC' + cluster_id].add(sample_id)

    # Read candidate genes per patient
    candidate_genes = {}  # sample_id -> { (gene_name, gene_type, diagnosis, pheno_match) }
    with open(candidate_genes_list) as chandle:
        header = None
        for cline in chandle:
            if cline.startswith('#'):
                continue
            if not header:
                header = cline.strip().split('\t')
                continue
            v = cline[:-1].split('\t')
            if v[header.index('patient_id')] not in candidate_genes:
                candidate_genes[v[header.index('patient_id')]] = set()
            candidate_genes[v[header.index('patient_id')]].add((
                v[header.index('gene_name')],
                v[header.index('gene_type')] if 'gene_type' in header else '',
                v[header.index('diagnosis')] if 'diagnosis' in header else '',
                v[header.index('pheno_match')] if 'pheno_match' in header else '',
            ))

    outhandle = open(cluster_details, 'w')
    outhandle.write(
        '# Clusters created using agglomerative clustering on normalized all-against-all pairwise HPO Phrank similarity\n' +
        '\t'.join(['cluster_id', 'patient_id', 'gene_name', 'gene_type', 'diagnosis', 'pheno_match'])
    )
    for cluster_id in sorted(cluster_to_patients.keys()):
        current_genes = set()
        for sample_id in cluster_to_patients[cluster_id]:
            for gene in candidate_genes.get(sample_id, {}):
                current_genes.add(gene[0])
        outhandle.write(('# Cluster {cluster_id} has {patient_cnt} total patients, {patient_candidates} ' +
                         'patients with candidate genes, and {gene_count} unique genes.\n').format(
            cluster_id=cluster_id,
            patient_cnt=len(cluster_to_patients[cluster_id]),
            patient_candidates=len([pt for pt in cluster_to_patients[cluster_id]
                                    if pt in candidate_genes and len(candidate_genes[pt]) > 0]),
            gene_count=len(current_genes)
        ))
        for sample_id in cluster_to_patients[cluster_id]:
            for gene in candidate_genes.get(sample_id, {}):
                outhandle.write('\t'.join([
                    cluster_id,
                    sample_id,
                    gene[0],  # gene_name
                    gene[1],  # gene_type
                    gene[2],  # diagnosis
                    gene[3]  # phenotypic match
                ])+'\n')
    outhandle.close()


########################################################################################################

def get_query_genes(cluster_details=cfg.script_directory+'/test/phrank/test_cluster_details.tsv'):
    """
    :param cluster_details: tab-delimited file produced by "create_gene_groups" function
    :return: set of genes to be provided to g:Profiler API call
    """

    query = {}

    with open(cluster_details) as chandle:
        header = None
        for cline in chandle:
            if cline.startswith('#'):
                continue
            if not header:
                header = cline.strip().split('\t')
                continue
            v = cline[:-1].split('\t')

            # add a "c" prefix to everything to force cluster_ids to be non-numerical
            if 'c' + v[header.index('cluster_id')] not in query:
                query['c' + v[header.index('cluster_id')]] = set()
            query['c' + v[header.index('cluster_id')]].add(v[header.index('gene_name')])

    return {cluster_id: sorted(list(genes)) for cluster_id, genes in query.items()}


########################################################################################################

def candidate_gene_info(cluster_details=cfg.script_directory+'/test/phrank/test_cluster_details.tsv'):
    """
    :param cluster_details: tab-delimited file produced by "create_gene_groups" function
    :return: string indicating which genes from which patients were implicated (for easier output parsing)
    """

    compelling_genes = {}
    with open(cluster_details) as chandle:
        header = None
        for cline in chandle:
            if cline.startswith('#'):
                continue
            if not header:
                header = cline.strip().split('\t')
                continue
            v = cline[:-1].split('\t')

            parenthetical = v[header.index('gene_type')]  # we'll return the "type" of gene if there is one specified

            if True in [value in v[header.index('diagnosis')].lower() for value in ['yes', 'y', 'diagnosis', 'true', '']]:
                parenthetical += ('|' if len(parenthetical) > 0 else '')+'DIAG'
            elif float(v[header.index('pheno_match')]) > 5:
                parenthetical +=('|' if len(parenthetical) > 0 else '') + str(round(float(v[header.index('pheno_match')]), 2))

            if 'c' + v[header.index('cluster_id')] not in compelling_genes:
                compelling_genes['c' + v[header.index('cluster_id')]] = {}
            compelling_genes['c' + v[header.index('cluster_id')]][v[header.index('gene_name')]] = v[header.index(
                'gene_name')] + '(' +parenthetical+')'

    return compelling_genes


########################################################################################################

def run_gene_set_enrichment(cluster_details=cfg.script_directory+'/test/phrank/test_cluster_details.tsv',
                            cluster_enrichments=cfg.script_directory+'/test/phrank/test_cluster_enrichments.tsv'):
    """
    :param cluster_details: tab-delimited file produced by "create_gene_groups" function
    :param cluster_enrichments: tab-delimited results from running g:Profiler's gene set enrichment analysis
    :return: (None) tab-delimited results from running g:Profiler's gene set enrichment analysis
    """

    # Read in the set of query genes and corresponding details
    query_genes = get_query_genes(cluster_details)
    compelling_genes = candidate_gene_info(cluster_details)

    # Send the request to g:Profiler API
    r = requests.post(
        url='https://biit.cs.ut.ee/gprofiler/api/gost/profile/',
        json={
            'user_threshold': 0.99,
            'organism': 'hsapiens',
            'sources': ["KEGG", "REAC"],
            'domain_scope': 'annotated',
            'query': query_genes
        })

    # Get Ensembl mapping (to print out gene names rather than Ensembl IDs)
    ensg_mapping = {}
    for query in query_genes.keys():
        if query in r.json()['meta']['genes_metadata']['query']:
            for gene_name, ensg_list in r.json()['meta']['genes_metadata']['query'][query]['mapping'].items():
                for ensg in ensg_list:
                    ensg_mapping[ensg] = gene_name

    # Print out enrichment results
    outhandle = open(cluster_enrichments, 'w')
    outhandle.write('\t'.join([
        'pathway_id', 'pathway_name', 'pathway_size', 'cluster_id', 'query_size', 'p_value', 'core_genes'
    ]) + '\n')

    for result in r.json()['result']:

        if result['native'] in ['REAC:0000000', 'KEGG:00000']:  # skip root pathways
            continue

        if len([i for i, v in enumerate(result['intersections']) if len(v) > 0]) > 1:  # for each enriched pathway

            outhandle.write('\t'.join(map(str, [
                result['native'],
                result['name'],
                result['term_size'],
                result['query'],
                result['query_size'],
                result['p_value'],
                ','.join([
                    compelling_genes.get(result['query'], {}).get(
                        ensg_mapping[r.json()['meta']['genes_metadata']['query'][result['query']]['ensgs'][i]],
                        ensg_mapping[r.json()['meta']['genes_metadata']['query'][result['query']]['ensgs'][i]]) for i, v
                    in enumerate(result['intersections']) if len(v) > 0
                ])
            ])) + '\n')
    outhandle.close()


########################################################################################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Redo pathway analysis')
    parser.add_argument('--pairwise_similarity', dest='pairwise_similarity', action='store_true', default=False,
                        help='Compute all-against-all patient phenotypic similarities from list of per-patient ' +
                             'HPO terms.')
    parser.add_argument('--cluster_patients', dest='cluster_patients', action='store_true', default=False,
                        help='Create patient clusters using similarity scores computed with `--pairwise_similarity`.')
    parser.add_argument('--query_genes', dest='query_genes', action='store_true', default=False,
                        help='Create sets of query genes for each patient cluster generated with `--cluster_patients`.')
    parser.add_argument('--run_gsea', dest='run_gsea', action='store_true', default=False,
                        help='Run gene set enrichment analysis on query genes per patient cluster and write to file.')
    parser.add_argument('--hpo_file', type=str, default=None,
                        help='Full path to a tab-delimited file containing HPO terms per patient')
    parser.add_argument('--similarity_file', type=str, default=None,
                        help='Full path to a tab-delimited file containing all-against-all patient phenotypic ' +
                             'similarity scores')
    parser.add_argument('--cluster_assignments_file', type=str, default=None,
                        help='Full path to a tab-delimited file containing cluster IDs -> patients in cluster')
    parser.add_argument('--candidate_genes_file', type=str, default=None,
                        help='Full path to a tab-delimited file containing candidate genes per patient')
    parser.add_argument('--query_genes_file', type=str, default=None,
                        help='Full path to a tab-delimited file containing query genes for each patient cluster')
    parser.add_argument('--geneset_enrichments_file', type=str, default=None,
                        help='Full path to a tab-delimited file containing enriched gene sets found using g:Profiler')
    args = parser.parse_args()

    if args.pairwise_similarity:
        sys.stderr.write('Computing all-against-all pairwise patient similarities...\n')

        if not args.hpo_file or not args.similarity_file:
            sys.stderr.write('\n'.join([
                '[ERROR] No file specifying per-patient HPO terms OR output file with similarities. Please rerun:',
                'python pathway_analysis.py \\',
                '    --pairwise_similarity \\',
                '    --hpo_file=<full/path/to/patient_hpo_terms_file.tsv> \\',
                '    --similarity_file=<full/path/to/output_similarities_file.tsv>'
            ])+'\n')
            sys.exit(1)

        pairwise_phrank_similarity(args.hpo_file, args.similarity_file)

    if args.cluster_patients:
        sys.stderr.write('Clustering patients using agglomerative complete-linkage clustering...\n')

        if not args.similarity_file or not args.cluster_assignments_file:
            sys.stderr.write('\n'.join([
                '[ERROR] No file containing patient pairwise similarities OR output file with cluster assignments. Please rerun:',
                'python pathway_analysis.py \\',
                '    --cluster_patients \\',
                '    --similarity_file=<full/path/to/output_similarities_file.tsv> \\',
                '    --cluster_assignments_file=<full/path/to/patient_to_cluster_assignments_file.tsv>'
            ])+'\n')
            sys.exit(1)

        create_patient_clusters(args.similarity_file, args.cluster_assignments_file)

    if args.query_genes:
        sys.stderr.write('Generating sets of query genes per patient cluster...\n')

        if not args.cluster_assignments_file or not args.candidate_genes_file or not args.query_genes_file:
            sys.stderr.write('\n'.join([
                '[ERROR] No file(s) containing patient cluster assignments, candidate genes per patient OR output file with query gene sets. Please rerun:',
                'python pathway_analysis.py \\',
                '    --query_genes \\',
                '    --cluster_assignments_file=<full/path/to/patient_to_cluster_assignments_file.tsv> \\',
                '    --candidate_genes_file=<full/path/to/patient_to_candidate_genes_file.tsv> \\',
                '    --query_genes_file=<full/path/to/query_genes_per_cluster.tsv>'
            ])+'\n')
            sys.exit(1)

        create_gene_groups(args.cluster_assignments_file, args.candidate_genes_file, args.query_genes_file)

    if args.run_gsea:
        sys.stderr.write('Running g:Profiler\'s gene set enrichment analysis...\n')

        if not args.query_genes_file or not args.geneset_enrichments_file:
            sys.stderr.write('\n'.join([
                '[ERROR] No file(s) containing query gene sets per cluster OR output file with enrichment results. Please rerun:',
                'python pathway_analysis.py \\',
                '    --run_gsea \\',
                '    --query_genes_file=<full/path/to/query_genes_per_cluster.tsv> \\',
                '    --geneset_enrichments_file=<full/path/to/output/gene_set_enrichment_results.tsv'
            ])+'\n')
            sys.exit(1)

        run_gene_set_enrichment(args.query_genes_file, args.geneset_enrichments_file)
