to_set_py <- NULL
to_graph_py <- NULL
diamond_core <- NULL


.onLoad <- function(libname, pkgname){
  path <- system.file(package = "MODifieRDev")
  diamond_module <- reticulate::import_from_path(module = "DIAMOnD_MODifieR_rt", path = path)
  
  to_set_py <<- diamond_module$to_set
  to_graph_py <<- diamond_module$to_graph
  diamond_core <<- diamond_module$DIAMOnD
  
  utils::globalVariables(c("clique", "i", "j"))
  
}