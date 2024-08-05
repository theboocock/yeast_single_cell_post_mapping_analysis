# Processed data files for Boocock et al. 2024

## Description of the data and file structure

This folder contains all the files necessary to recreate the figures and analyses described and presented in [https://elifesciences.org/reviewed-preprints/95566](https://elifesciences.org/reviewed-preprints/95566).

For compatibility with the scripts on the associated github [https://github.com/theboocock/yeast_single_cell_post_mapping_analysis](https://github.com/theboocock/yeast_single_cell_post_mapping_analysis). When expanded the files should be extracted to a sub-directory name `data` within the directory containing the script files from github. The files described below are used by the code on github to reproduce the figures and tables presented in our manuscript on single-cell eQTL mapping in yeast.

* The file `data/bloom_causal_variant_mapping.csv` contains a list of causal genes for QTL as described in previous work from the same crosses we used for single-cell eQTL mapping [https://elifesciences.org/articles/49212](https://elifesciences.org/articles/49212). Each row contains the causal gene, associated statistics, and other genomic annotations.
* The file `data/bulkEQTL_cis_only_test.RDS` is an RDS object containing the result of performing a cis only test using a count-based model with the data described a 2018 paper from our lab that measured gene expression with RNA-seq [https://elifesciences.org/articles/35471](https://elifesciences.org/articles/35471).
* The file `data/cell_cycle_feb02162022.tsv` contains cell-cycle assignments for all cells used. Each row contains the cell-cycle assignment for cells analyzed in the experiment. Columns include the data set, cell barcode, cell-cycle assignment, and Seurat-based cluster assignments prior to manual cell-cycle assignment.
* The file `data/haploid_genes.csv` contains a list of haploid-specific genes in yeast curated from the literature.
* The file `data/mating1.csv` contains the results of one replicate of our mating assay comparing the 82R variant of GPA1 to the 82W variant.
* The file `data/mating2.csv` contains the second replicate of the experiment described above.
* The file `data/qtl_bloom.csv`  contains a list of causal genes for QTL as described in previous work from the same crosses we used for single-cell eQTL mapping. Each row contains a QTL, associated statitics, and other genomic annotations.
* The file `data/saccharomyces_cerevisiae.gff` is the genomic annotation file downloaded from SGD [https://www.yeastgenome.org/](https://www.yeastgenome.org/) on Jan 17, 2017 with the version R64-2-1.
* The file `data/single_cell_runs.csv` contains a list of all the single cell data sets used in our analyses. Each row contains a dataset name, the type of run, and whether the run passed an initial QC step. Note: runs with a N for the good were not analysed in our manuscript.
* The file `data/cc/botstein_genes.csv` contains a list of the 799 cell-cycle genes described in [https://pubmed.ncbi.nlm.nih.gov/9843569/](https://pubmed.ncbi.nlm.nih.gov/9843569/). Each row contains the name of the systematic name of the gene.
* The file `data/cc/ccdf.20210511.RDS` is an RDS object containing a data frame of the 22 cell-cycle genes we curated.
* The file `data/cc/cell_cycle.csv` contains a list of the 799 cell-cycle genes described in [https://pubmed.ncbi.nlm.nih.gov/9843569/](https://pubmed.ncbi.nlm.nih.gov/9843569/). Each row contains the systematic name of the gene, a score, a phase, and the assigned peak cell-cycle for the expression of the gene. This information was helpful when annotating the cell-cycle stage of the clusters identified using Seurat.
* The file `data/images/one_pot.png` contains the png image used in the one_pot diagram shown in figure S1.
* The file `data/joseph/1012.gds` contains the SNPrelate object with the genotype matrices of the the 1,011 strains from Peter et al. for all variants with a frequency of >5%. The 1,011 dataset is described in more detail here [https://pubmed.ncbi.nlm.nih.gov/29643504/](https://pubmed.ncbi.nlm.nih.gov/29643504/).
* The file `data/joseph/full_annotations.xlsx` and in particular the sheet `peter2018` contains the annotations for the 1,011 strains from Peter et al. . This data was taken from the supplementary information of [https://pubmed.ncbi.nlm.nih.gov/29643504/](https://pubmed.ncbi.nlm.nih.gov/29643504/). Any cell with #N/A as its value indicates that data was not available.
* The file `data/platereader/galmutantsglufeb2019conditions.csv` is an annotation file used with the script `load_gpa1_phenotype_and_mating_data.R` to analyze plate reader data.
* The folder `data/plate_reader/gpa1_validation/*` contains various files representing the plate reader data for our growth experiments comparing the growth of strains with the 82R variant of GPA1 to strains with the 82W variant. The script `load_gpa1_phenotype_and_mating_data.R` contains the logic to process these files. 
* The file `data/provean/1002codingEffects.RDS` contains the PROVEAN scores for all coding variants identified in the 1,011 yeast genomes study.
* The folder `data/ref/yeast/` contains the STAR reference file used with the software cellranger to analyse all single-cell data.
* The file `data/rr/extracted_average_phenotypes.RData` is an RData object that contains the growth measurements for all haploid segregants from our large-scale growth trait mapping study that explored the same crosses we used for our single-cell eQTL mapping.
* The file `data/rr/seg.recoded.RData` is an RData object with the genotypes for all haploid segregants from our large-scale growth trait mapping study that explored the same crosses we used for our single-cell eQTL mapping.
* The folder `single_cell` contains the single-cell used for validating the impact of the 82R variant of GPA1 on gene expression and the distribution of cell-cycle stages. The analysis of this data was done with cellranger.
* The folder `vcf/all_no_filt.txt` contains the unfiltered variants and their effects predicted with SNPEff for the haploid parental strains we sequenced in our large-scale trait mapping study. This unfiltered dataset did not exclude indels so to make it as useful as possible for annotated each trans-eQTL hotspot. Each row contains a different variant and associated metadata.
* The folder `data/out/*` contains subdirectories for each single-cell dataset. These folders contain the data generated from the HMM and mapping analysis code. These files in these directories are necessary to reproduce the cell-cycle annotation analysis, and generate the figures and tables presented in our manuscript.
* The folder `data/out/combined/*` contains the merged analysis of the crosses used for single-cell eQTL mapping. Each folder is named based on the cross identifier. The files in these directories are necessary to generate the figures and tables of our manuscript.
* The folder `data/out/cell_cycle/*` contains the cell-cycle annotations and associated RDS objects for each single-cell dataset.

Other folders and files in these the `data/out/` may not be used by the final version of the analysis code, but we included them for record keeping purposes as they contain information that may be useful for researchers following up on our work.

## Sharing/Access information

Much of this data is derived from the initial HMM and mapping analysis. The data dryad link and associated github is provided below.

* [https://doi.org/10.5061/dryad.xgxd254qb](https://doi.org/10.5061/dryad.xgxd254qb)
* [https://github.com/joshsbloom/single_cell_eQTL/tree/master/yeast/code](https://github.com/joshsbloom/single_cell_eQTL/tree/master/yeast/code)

The annotation data was obtained from various sources that are cited throughout the manuscript linked below.

* [https://elifesciences.org/reviewed-preprints/95566](https://elifesciences.org/reviewed-preprints/95566)

## Code/Software

Analysis code for the post mapping analyses can provided be found at [https://github.com/theboocock/yeast_single_cell_post_mapping_analysis](https://github.com/theboocock/yeast_single_cell_post_mapping_analysis)

