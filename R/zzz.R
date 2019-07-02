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