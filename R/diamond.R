#' DIAMOnD module
#' 
#' A seed gene based algorithm to identify disease module from differentially expressed genes
#' 
#' @inheritParams clique_sum
#' @param n_output_genes maximum number of genes to be included in the final module 
#' @param seed_weight Numeric additional parameter to assign weight for the seed genes 
#' @param include_seed Logical TRUE/FALSE for inclusion of seed genes in the output module 
#' @details
#' A slightly modified version of the original DIAMOnD python script is called from within R.
#' The only change to the orginal algorithm is the option to include the seed genes
#' to the module. There are also function to add or remove the seed genes from the output object, namely:
#' \code{\link{diamond_add_seed_genes}} and \code{\link{diamond_remove_seed_genes}}
#' For a detailed description of how the algorithm works, please see the paper referenced below.
#' @return 
#' diamond returns an object of class "MODifieR_module" with subclass "DIAMOnD". 
#' This object is a named list containing the following components:
#' \item{module_genes}{A character vector containing the genes in the final module}
#' \item{seed_genes}{Character vector containing genes that have been used as seed genes in the algorithm}
#' \item{ignored_genes}{Potential seed genes that are not in the PPi network}
#' \item{added_genes}{A table containing information on all added genes. First column is the name of the gene,
#' the second column is the degree of the node (gene). The third column is the number of +1 neighbors 
#' and the fourth column is the p-value.}
#' \item{settings}{A named list containing the parameters used in generating the object}
#' @references 
#' \cite{Ghiassian, S. D., Menche, J., & Barabási, A. L. (2015). 
#' A DIseAse MOdule Detection (DIAMOnD) Algorithm Derived from a Systematic Analysis of 
#' Connectivity Patterns of Disease Proteins in the Human Interactome. PLoS Computational 
#' Biology, 11(4), 1–21. \url{https://doi.org/10.1371/journal.pcbi.1004120}}
#' @export
diamond <- function(MODifieR_input, ppi_network, deg_cutoff = 0.05, n_output_genes = 200, seed_weight = 10,
                               include_seed = FALSE, dataset_name = NULL){
  # Retrieve settings
  settings <- do.call(what = "settings_function", as.list(stackoverflow::match.call.defaults()[-1]))
  
  if (!is.null(dataset_name)){
    settings$MODifieR_input <- dataset_name
  }
  
  #Get the input genes
  diamond_genes <- MODifieR_input$diff_genes
  diamond_genes <- diamond_genes[diamond_genes$pvalue < deg_cutoff, ]
  diamond_genes <- unique(stats::na.omit(diamond_genes$gene)) 
  #Convert to python objects
  diamond_set <- to_set_py(diamond_genes)
  ppi_network <- as.matrix(ppi_network)
  storage.mode(ppi_network) <- "character"
  ppi_graph <- to_graph_py(as.matrix(ppi_network), nx$Graph())
  #Run python scri[t]
  raw_module_py <- diamond_core(ppi_graph, diamond_set, as.integer(n_output_genes), as.integer(seed_weight))
  #Extract data
  raw_module <- raw_module_py[[1]]
  seed_genes <- raw_module_py[[2]] 
  ignored_genes <- raw_module_py[[3]]
  #Prepare the added genes info table
  added_genes <- data.frame(cbind(sapply(X = 1:4, FUN = function(i){sapply(raw_module, function(x)x[[i]])})), stringsAsFactors = F)
  added_genes[,2:4] <- sapply(added_genes[,2:4], FUN = as.numeric)
  colnames(added_genes) <-   c("Gene","Degree", "Connectivity", "p-value")
  
  module_genes <- sapply(raw_module, function(x)x[[1]])
  if (include_seed){
    module_genes <- c(module_genes, seed_genes)
  }
  # Build new MODifieR object
  new_diamond_module <- list("module_genes" = module_genes,
                             "seed_genes" =  seed_genes,
                             "ignored_genes" = ignored_genes,
                             "added_genes" = added_genes,
                             "settings" = settings)
  
  class(new_diamond_module) <- c("MODifieR_module", "DIAMOnD")
  
  return(new_diamond_module)
}

#' Remove seed genes from a DIAMOnD MODifieR_module
#' @param diamond_module A \code{MODifieR_input} object created by \code{\link{diamond}}
#' @details 
#' Removes seed genes from a DIAMOnD module
#' @return
#' An object of class "MODifieR_module" with subclass "DIAMOnD"
#' @seealso 
#' \code{\link{diamond}}
#' @export
diamond_remove_seed_genes <- function(diamond_module){
  diamond_module$module_genes <- setdiff(diamond_module$module_genes, diamond_module$seed_genes)
  diamond_module$settings$include_seed <- F
  return(diamond_module)
}
#' Add seed genes from a DIAMOnD MODifieR_module
#' @inheritParams diamond_remove_seed_genes
#' @details 
#' Adds seed genes from a DIAMOnD module
#' @return
#' An object of class "MODifieR_module" with subclass "DIAMOnD"
#' @seealso 
#' \code{\link{diamond}}
#' @export
#' @export
diamond_add_seed_genes <- function(diamond_module){
  diamond_module$module_genes <- unique(c(diamond_module$module_genes, diamond_module$seed_genes))
  diamond_module$settings$include_seed <- T
  return(diamond_module)
}