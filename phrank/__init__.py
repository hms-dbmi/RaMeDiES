"""
Downloaded from Bitbucket:
https://bitbucket.org/bejerano/phrank/raw/7e0adbbe07cb12589ddba730565d283f6f06507e/phrank/__init__.py

Phrank publication: https://pubmed.ncbi.nlm.nih.gov/29997393/
"""

from collections import defaultdict
import math
from .utils import load_maps, load_term_hpo, closure, load_disease_gene, compute_gene_disease_pheno_map


class Phrank:
    @staticmethod
    def compute_information_content(annotations_map, child_to_parent_map):
        information_content, marginal_information_content = {}, {}
        annotatedGeneCt = 0
        associated_phenos = defaultdict(set)
        for gene, phenos in annotations_map.items():
            annotatedGeneCt += 1
            all_ancestors = closure(phenos, child_to_parent_map)

            # for each ancestor increment the count since this pheno is now associated with the specified gene
            for pheno in all_ancestors:
                associated_phenos[pheno].add(gene)

        phenos = associated_phenos.keys()
        for pheno in phenos:
            information_content[pheno] = -math.log(1.0 * len(associated_phenos[pheno]) / annotatedGeneCt, 2) if len(
                associated_phenos[pheno]) else 0

        for pheno in phenos:
            parent_phenos = child_to_parent_map[pheno]
            parent_entropy = 0
            if len(parent_phenos) == 1:
                parent_entropy = information_content[parent_phenos[0]]
            elif len(parent_phenos) > 1:
                list_of_phenosets = [associated_phenos[parent] for parent in parent_phenos]
                parent_set = set([])
                for phenoset in list_of_phenosets:
                    parent_set = parent_set | phenoset if parent_set else phenoset
                parent_entropy = -math.log(1.0 * len(parent_set) / annotatedGeneCt, 2) if len(parent_set) else 0
            marginal_information_content[pheno] = information_content[pheno] - parent_entropy
        return information_content, marginal_information_content

    def __init__(self, dagfile, diseaseannotationsfile=None, diseasegenefile=None, geneannotationsfile=None):
        """Initialize Phrank object with the disease annotations file or gene annotations file"""
        self._child_to_parent, self._parent_to_children = load_maps(dagfile)
        if diseaseannotationsfile and diseasegenefile:
            self._disease_pheno_map = load_term_hpo(diseaseannotationsfile)
            self._disease_gene_map = load_disease_gene(diseasegenefile)
            self._gene_pheno_map = compute_gene_disease_pheno_map(self._disease_gene_map, self._disease_pheno_map)
            self._IC, self._marginal_IC = Phrank.compute_information_content(self._gene_pheno_map,
                                                                             self._child_to_parent)
            self._gene_and_disease = True
        elif geneannotationsfile:
            self._gene_pheno_map = load_term_hpo(geneannotationsfile)
            self._IC, self._marginal_IC = Phrank.compute_information_content(self._gene_pheno_map,
                                                                             self._child_to_parent)
            self._gene_and_disease = False

    def get_term_ic(self, phenotype):
        return self._IC.get(phenotype, 0)

    def get_causal_rank(self, scores, causal_item):
        rank = 1
        for score in scores:
            if score[1] in causal_item:
                return score[1], score[0], rank
            else:
                rank = rank + 1
        return causal_item, 0, len(scores) + 1

    def rank_diseases(self, patient_genes, patient_phenotypes, baseline=False):
        """Compute the Phrank score for each disease matching the patient phenotypes"""
        disease_scores = []
        for disease in self._disease_pheno_map:
            if self._disease_gene_map[disease] & patient_genes:
                disease_phenos = self._disease_pheno_map.get(disease, set([]))
                score = self.compute_phenotype_match(patient_phenotypes, disease_phenos)
                score = self.compute_phenotype_match(patient_phenotypes,
                                                     disease_phenos) if not baseline else self.compute_baseline_match(
                    patient_phenotypes, disease_phenos)
                disease_scores.append((score, disease))
        disease_scores.sort(reverse=True)
        return disease_scores

    def rank_genes(self, patient_genes, patient_phenotypes, normalized=False, baseline=False):
        gene_scores = []
        if self._gene_and_disease:
            gene_scores = self.rank_genes_using_disease(patient_genes, patient_phenotypes, normalized=normalized,
                                                        baseline=baseline)
        else:
            gene_scores = self.rank_genes_directly(patient_genes, patient_phenotypes, normalized=normalized,
                                                   baseline=baseline)
        return gene_scores

    def rank_genes_directly(self, patient_genes, patient_phenotypes, normalized=False, baseline=False):
        gene_scores = []
        for gene in patient_genes:
            gene_phenos = self._gene_pheno_map[gene]
            score = self.compute_phenotype_match(patient_phenotypes,
                                                 gene_phenos) if not baseline else self.compute_baseline_match(
                patient_phenotypes, gene_phenos)
            if normalized:
                max_gene_score = self.compute_maximal_match(gene_phenos)
                score = 1.0 * score / max_gene_score
            gene_scores.append((score, gene))
        gene_scores.sort(reverse=True)
        return gene_scores

    def rank_genes_using_disease(self, patient_genes, patient_phenotypes, normalized=False, baseline=False):
        """Compute the Phrank score for each gene matching the patient phenotypes"""
        genedisease_scores = defaultdict(list)
        for disease in self._disease_pheno_map:
            if self._disease_gene_map[disease] & patient_genes:
                disease_phenos = self._disease_pheno_map.get(disease, set([]))
                score = self.compute_phenotype_match(patient_phenotypes,
                                                     disease_phenos) if not baseline else self.compute_baseline_match(
                    patient_phenotypes, disease_phenos)
                if normalized:
                    max_disease_score = self.compute_maximal_disease_match(disease)
                    score = 1.0 * score / max_disease_score
                for gene in self._disease_gene_map[disease] & patient_genes:
                    genedisease_scores[gene].append(score)

        gene_scores = []
        for gene in genedisease_scores:
            score = max(genedisease_scores[gene])
            gene_scores.append((score, gene))
        gene_scores.sort(reverse=True)
        return gene_scores

    def compute_maximal_match(self, phenotypes):
        """
        input: a list of phenotypes [HPO:XXX, HPOYYY]
        output: score (Float) - the maximal score this set of phenotypes can receive

        The maximum score achievable for a set of phenotypes is if every phenotype
        is found in the comparing set. This number can be used to normalize the phenotype
        score so that query with a large # of phenotypes do not skew the results
        """
        return self.compute_phenotype_match(phenotypes, phenotypes)

    def compute_maximal_disease_match(self, disease):
        disease_phenos = self._disease_pheno_map.get(disease, set([]))
        return self.compute_phenotype_match(disease_phenos, disease_phenos)

    def compute_phenotype_match(self, patient_phenotypes, query_phenotypes):
        """
        input: patient_phenotypes, query phenotypes - two lists of phenotypes
        output: score - Phrank score measuring the similarity between the two sets
        """
        # change_to_primary, patient_genes, disease_gene_map, disease_pheno_map, child_to_parent, disease_marginal_content)
        all_patient_phenotypes = closure(patient_phenotypes, self._child_to_parent)
        all_query_phenotypes = closure(query_phenotypes, self._child_to_parent)
        similarity_score = 0
        for phenotype in all_patient_phenotypes & all_query_phenotypes:
            similarity_score += self._marginal_IC.get(phenotype, 0)
        return similarity_score

    def compute_baseline_match(self, patient_phenotypes, query_phenotypes):
        all_patient_phenotypes = closure(patient_phenotypes, self._child_to_parent)
        all_query_phenotypes = closure(query_phenotypes, self._child_to_parent)
        return len(all_patient_phenotypes & all_query_phenotypes)