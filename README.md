# Scripts for "Single-cell eQTL mapping in yeast reveals a tradeoff between growth and reproduction" 
## This repository also contains the code and data to recreate the figures, tables, performs trans-eQTL hotspot analysis, cell-cycle stage assignment, raw data processing, and additional links to generated data.

Sequencing data for each single-cell experiment is available on the SRA [PRJNA1049497](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1049497).

see [load_all_data.R](load_all_data.R) for the main post analysis script that generates the figures and tables, and performs the trans-eQTL hotspot analysis.

To run these scripts one needs to first download the processed data from Dryad and extract it into the data directory. the archive on dryad was created by splitting a large archive to recreate the original archive simply concatenate them together using the command `cat additional_data.tar.p* > additional_data.tar.gz`. The data can be downloaded [here](https://datadryad.org/stash/share/vHoXP5Ph78EPD3C2yGvnuV4F1DukxYyS9dUqa75_Wkg). 

A detailed description of the structure of the files in the archive can be found at the Dryad link. 

 
-----

see [cell_cycle_annotation.R](cell_cycle_annotation/cell_cycle_annotation.R) for the cell-cycle stage assignment script.

-----------------

see [numbers_final](https://github.com/theboocock/yeast_single_cell_post_analysis/tree/main/numbers_final) for scripts that generate all the numbers found in the maintext.

----

see [raw_data](raw_data_processing/extract_parents_and_vatrix_hoff.sh) for the raw data processing script.

see [pre_processing_repository](https://github.com/joshsbloom/single_cell_eQTL/tree/master/yeast/code) to access the code used to run the HMM and eQTL mapping analysis.

------------------

see [R_dependencies](R_dependencies.yaml) for the R version and list of packages used for this analysis.
