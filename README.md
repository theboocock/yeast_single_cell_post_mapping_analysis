# Scripts for "Single-cell eQTL mapping in yeast reveals a tradeoff between growth and reproduction" 
## This repository also contains the code that recreates the figures, tables, performs trans-eQTL hotspot analysis, cell-cycle stage assignment, raw data processing, and additional links to generated data.

Raw data folder for the post analysis scripts are located [here](https://drive.google.com/drive/folders/1SAUYxO7EhUq-FQLzrc__Lm_0dVF06oIj?usp=drive_link)

sequencing data for each single-cell experiment is available on the SRA [PRJNA](http://ncbi.com)

see [load_all_data](load_all_data.R) for the main post analysis script that generates the figures and tables, and performs the trans-eQTL hotspot analysis.

-----

see [cell_cycle_annotation](cell_cycle_annotation/cell_cycle_annotation.R) for the cell-cycle stage assignment script

-----------------

see [raw_data](raw_data_processing/extract_parents_and_vatrix_hoff.sh) for the raw data processing script

------------------

see [drivelink](drivelink.com) for the genotype and gene expression matrices, and the cell meta-data for each experiment in a simple format.
