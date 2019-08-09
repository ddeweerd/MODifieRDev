to_set_py <- NULL
to_graph_py <- NULL
diamond_core <- NULL
build_clique_db_py <- NULL
nx <- NULL
np <- NULL
scipy <- NULL
sqlite <- NULL

unevaluated_args <- c("MODifieR_input", "MODifieR_module1","MODifieR_module2",
                      "ppi_network", "module_list", "expression_matrix", 
                      "annotation_table", "count_matrix")

unexported_settings <- c("group1_indici", "group2_indici")

exported_rownames <- c("annotated_exprs_matrix", "expression_matrix", "module_p_values", 
                       "edgeR_deg_table", "count_matrix")

exported_names <- c("color_vector", "gene_by_method")

remove_headers <- c("module_genes", "ignored_genes", "jaccard_table", "module_colors", "seed_genes",
                    "group1_modules", "group2_modules")

unit_intervals <- c("clique_significance", "deg_cutoff", "frequency_cutoff", "fraction_of_interactions",
                    "vwp", "pval_cutoff", "maxPOutliers")

integer_variables <- c("min_clique_size", "min_deg_in_clique", "n_cores", "n_iterations", "iteration",
                       "n_output_genes", "seed_weight", "hierarchy", "permutations", "repeats",
                       "minModuleSize", "deepSplit", "nrep", "nrep2")

boolean_variables <- c("multiple_cores", "include_seed", "haircut", "fluff", "loops", "pamRespectsDendro",
                       "saveTOMs", "numericLabels", "expression", "differential_expression", "filter_expression",
                       "use_adjusted", "normalize_quantiles")

numeric_variables <- c("module_score", "fdt", "specificTheta", "conservedTheta")


.onLoad <- function(libname, pkgname){
  nx <<- reticulate::import(module = "networkx", as = "nx", delay_load = T)
  np <<- reticulate::import(module = "numpy", as = "np", delay_load = T)
  scipy <<- reticulate::import(module = "scipy", delay_load = T)
  sqlite <<- reticulate::import(module = "sqlite3", delay_load = T)
  
  
  path <- system.file(package = "MODifieRDev")
  diamond_module <- reticulate::import_from_path(module = "DIAMOnD_MODifieR_rt", path = path)
  
  path <- system.file(package = "MODifieRDev")
  clique_db_module <- reticulate::import_from_path(module = "db_cliques", path = path)
  
  to_set_py <<- diamond_module$to_set
  to_graph_py <<- diamond_module$to_graph
  diamond_core <<- diamond_module$DIAMOnD
  build_clique_db_py <<- clique_db_module$build_clique_db
  
  utils::globalVariables(c("clique", "i", "j"))
  
}