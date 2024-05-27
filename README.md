# Scripts for "Single-cell eQTL mapping in yeast reveals a tradeoff between growth and reproduction" 
## This repository also contains the code and data to recreate the figures, tables, performs trans-eQTL hotspot analysis, cell-cycle stage assignment, raw data processing, and additional links to generated data.

Raw data folder for the post analysis scripts are located [here](https://drive.google.com/drive/folders/1SAUYxO7EhUq-FQLzrc__Lm_0dVF06oIj?usp=drive_link).

sequencing data for each single-cell experiment is available on the SRA [PRJNA1049497](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1049497).

see [load_all_data](load_all_data.R) for the main post analysis script that generates the figures and tables, and performs the trans-eQTL hotspot analysis.

-----

see [cell_cycle_annotation](cell_cycle_annotation/cell_cycle_annotation.R) for the cell-cycle stage assignment script.

-----------------

see [raw_data](raw_data_processing/extract_parents_and_vatrix_hoff.sh) for the raw data processing script.

see [pre_processing_repository](https://github.com/joshsbloom/single_cell_eQTL/tree/master/yeast/code) to access the code used to run the HMM and eQTL mapping analysis, from which one can obtain the genotype and gene expression matrices.

------------------

see [R_dependencies](R_dependencies.yaml) for the R version and list of packages used for this analysis.
