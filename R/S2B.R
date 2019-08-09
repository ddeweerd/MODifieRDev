#' S2B
#' 
#' Combine two different disease modules into using Specific Betweenness (S2B)
#' @inheritParams clique_sum_exact
#' @inheritParams build_clique_db
#' @param MODifieR_module1 A MODifieR module inferred by a disease module inference method
#' @param MODifieR_module2 A MODifieR module inferred by a disease module inference method
#' @param nrep number of randomizations (shuffle edges maintaining node degree) 
#' to compute specificity score 1
#' @param nrep2 number of randomizations (shuffle seed identity) to compute specificity score 2
#' 
#' @details 
#' S2B prioritizes genes frequently and specifically present in shortest paths linking two disease modules. 
#' For a detailed description of the algorithm, see the paper referenced below
#' 
#' @return 
#' Returns an object of class "MODifieR_module" with subclass "S2B". 
#' This object is a named list containing the following components:
#' \item{module_genes}{A character vector containing the genes in the final module}
#' \item{input1_class}{Method that was used to create disease module 1 (inferred by its class)}
#' \item{input2_class}{Method that was used to create disease module 2 (inferred by its class)}
#' \item{settings}{A named list containing the parameters used in generating the object}
#' @references 
#' \cite{Garcia-Vaquero, M. L., Gama-Carvalho, M., Rivas, J. D. Las, & Pinto, F. R. (2018). 
#' Searching the overlap between network modules with specific betweeness (S2B) and its application 
#' to cross-disease analysis. Scientific Reports, 8(1), 1–10. 
#' \url{https://doi.org/10.1038/s41598-018-29990-7}}
#' 
#' 
#' @export
S2B <- function(MODifieR_module1, MODifieR_module2, ppi_network, nrep = 100, nrep2 = 100, n_cores = 4){
  
  # Retrieve settings
  evaluated_args <- c(as.list(environment()))
  settings <- as.list(stackoverflow::match.call.defaults()[-1])
  replace_args <- names(settings)[!names(settings) %in% unevaluated_args]
  for (argument in replace_args) {
    settings[[which(names(settings) == argument)]] <- evaluated_args[[which(names(evaluated_args) == 
                                                                              argument)]]
  }
  
  if(class(MODifieR_module1)[1] != "MODifieR_module"){
    stop("MODifieR_module1 is not of class MODifieR_module", call. = F)
  }
  
  if(class(MODifieR_module2)[1] != "MODifieR_module"){
    stop("MODifieR_module2 is not of class MODifieR_module", call. = F)
  }
  
  validate_ppi(ppi_network)
  
  validate_inputs(settings)
  
  ppi_graph <- simpmain(igraph::graph_from_data_frame(ppi_network))
 
  set1_index <- seedrows(ppi_graph, MODifieR_module1$module_genes)
  set2_index <- seedrows(ppi_graph, MODifieR_module2$module_genes)
  
  s2b_results <- S2B_core(ppi_graph, set1_index, set2_index, nrep, nrep2, n_cores = 10)
  
  s2b_cutoff <- s2bthreshold(s2b_results$s2btable$S2B)
  
  s2b_genes <- s2b_results$s2btable$id[s2b_results$s2btable$S2B > s2b_cutoff ]
  
  s2b_genes <- as.vector(as.character(s2b_genes))
  
  module_class1 <- class(MODifieR_module1)[2]
  module_class2 <- class(MODifieR_module1)[2]
  
  new_s2b_object <- list("module_genes" = s2b_genes,
                         "input1_class" = class(MODifieR_module1)[2],
                         "input2_class" = class(MODifieR_module2)[2],
                         "settings" = settings)
  
  class(new_s2b_object) <- c("MODifieR_module", "S2B")
  
  return(new_s2b_object)
  
}