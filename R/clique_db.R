#' Create a sqlite database containing all maximal cliques from a network
#'
#' @inheritParams clique_sum_permutation
#' @param ppi_network A network as a dataframe where the first 2 columns are the interactions
#' @param db_folder A directory where the database will be stored
#' @param db_name The name of the database. File suffix ".sqlite" will be appended
#' @details
#' Creates a SQLite database containing all maximal cliques from the network in the folder 
#' \code{db_folder} with filename \code{db_name.sqlite}. This database can be used as in 
#' input to \code{\link{clique_sum_exact}} and \code{\link{clique_sum_permutation}}
#' 
#' @export
build_clique_db <- function(ppi_network, db_folder, db_name){
  db_name <- paste0(db_name, ".sqlite")
  graphed_frame <- igraph::graph.data.frame(unique(ppi_network) , directed = FALSE)
  graphed_frame <- igraph::simplify(graphed_frame)
  vertex_dictionary <- cbind(as.vector((V(graphed_frame)) -1), as.vector(names(V(graphed_frame))))
  utils::write.table(x = vertex_dictionary, 
              paste0(db_folder, "vertex_dictionary.txt"), 
              quote = F, 
              row.names = F, 
              col.names = F)
  infer_cliques(graphed_frame = graphed_frame,
                db_folder = db_folder, 
                db_name = db_name, 
                temp_cliques = paste0(db_folder, "temp_cliques.txt"),
                vertex_dictionary_file = paste0(db_folder, "vertex_dictionary.txt"))

  file.remove(paste0(db_folder, "temp_cliques.txt"))
  file.remove(paste0(db_folder, "vertex_dictionary.txt"))
}

infer_cliques <- function(graphed_frame, db_folder, db_name, temp_cliques, vertex_dictionary_file){
  igraph::max_cliques(graph = graphed_frame, file = temp_cliques)
  build_clique_db_py(db_folder, db_name, temp_cliques, vertex_dictionary_file)
}

