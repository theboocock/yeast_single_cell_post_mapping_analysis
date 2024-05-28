source("pipeline/load_cross_objects.R")
## Load this sessionS
#save.image("/media/theboocock/scratch/Rsession_set2.RData")
source("r_final/make_figures_and_tables.R")

devtools::session_info() %>% yaml::write_yaml("R_dependencies.yaml")
