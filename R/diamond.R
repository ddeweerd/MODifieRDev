#' DIAMOnD module
#' 
#' A seed gene based algorithm to identify disease module from Differentially Expressed Genes
#' 
#' @inheritParams clique_sum
#' @param n_output_genes maximum number of genes to be included in the final module 
#' @param seed_weight Numeric additional parameter to assign weight for the seed genes 
#' @param include_seed Logical TRUE/FALSE for inclusion of seed genes in the output module 
#' @details
#' A slightly modified version of the original DIAMOnD python script is called from within R. the exchange of 
#' data is handled by temporary files. The only change to the orginal algorithm is the option to include the seed genes
#' to the module. There are also function to add or remove the seed genes from the output object, namely:
#' \code{\link{add_diamond_seed_genes}} and \code{\link{remove_diamond_seed_genes}}
#' 
#' the only change is the option
#' to include seed genes in the final module. The original python script is called from within R and  For a detailed description of how the algorithm works, please see
#' the paper referenced below.
#' @return 
#' diamond returns an object of class "MODifieR_module" with subclass "DIAMOnD". 
#' This object is a named list containing the following components:
#' \item{module_genes}{A character vector containing the genes in the final module}
#' \item{seed_genes}{Character vector containing genes that have been used as seed genes in the algorithm}
#' \item{ignored_genes}{Potential seed genes that have not been used to infer the module}
#' \item{settings}{A named list containing the parameters used in generating the object}
#' @references 
#' \cite{Ghiassian, S. D., Menche, J., & Barabási, A. L. (2015). 
#' A DIseAse MOdule Detection (DIAMOnD) Algorithm Derived from a Systematic Analysis of 
#' Connectivity Patterns of Disease Proteins in the Human Interactome. PLoS Computational 
#' Biology, 11(4), 1–21. https://doi.org/10.1371/journal.pcbi.1004120}
#' @export
diamond <- function(MODifieR_input, ppi_network, deg_cutoff = 0.05, n_output_genes = 200, seed_weight = 10,
                    include_seed = T, tempfile_genes = tempfile(), dataset_name = NULL){
  # Retrieve settings
  default_args <- formals()
  user_args <- as.list(match.call(expand.dots = T)[-1])
  settings <- c(user_args, default_args[!names(default_args) %in% names(user_args)])
  
  if (!is.null(dataset_name)){
    settings$MODifieR_input <- dataset_name
  }
  
  diamond_genes <- MODifieR_input$diff_genes
  
  diamond_genes <- diamond_genes[diamond_genes$pvalue < deg_cutoff, ]
  diamond_genes <- unique(na.omit(diamond_genes$ENTREZID)) 

  input_genes <- tempfile_genes
  input_ppi <- tempfile()
  write.table(x = diamond_genes, file = input_genes, sep = "\t" ,
              row.names = FALSE , quote = FALSE, col.names = F)
  
  write.table(x = ppi_network, file = input_ppi, sep = "\t" ,
              row.names = FALSE , quote = FALSE, col.names = F)
  
  # Sets python path
  python_path <- paste(system.file(package="MODifieRDev"), "DIAMOnD_MODifieR.py", sep="/")
  # Concetenates python call with CL arguments
  python_call <- paste("python", python_path, input_ppi, input_genes,
                       n_output_genes, seed_weight, sep = " ")
  # Gets output from python script as character vector
  raw_module <- system(command = python_call, intern = T)
  # Python scripts returns concatented ouput delimited by ",", split it
  split_module <- sapply(X = raw_module, FUN = split_values)
  # Include seed genes used?
  if (include_seed == T){
    module_genes <- c(split_module[[2]], split_module[[3]])
  } else{
    module_genes <- split_module[[3]]
  }
  # Build new MODifieR object
  new_diamond_module <- list("module_genes" = module_genes,
                             "seed_genes" =  split_module[[2]],
                             "ignored_genes" = split_module[[1]],
                             "added_genes" = as.data.frame(x = (split_module[3:6]),
                                                           row.names = 1:length(split_module[[3]]),
                                                           col.names =      c("Gene",
                                                                              "Degree",
                                                                              "Connectivity",
                                                                              "p-value"),
                                                           stringsAsFactors = F),
                             "settings" = settings)

  class(new_diamond_module) <- c("MODifieR_module", "DIAMOnD")

  return(new_diamond_module)
}
# Splits character vector on ","
split_values <- function(values){
  values <- gsub(pattern = " ", replacement = "", x = values)
  strsplit(x = values, split = ",")
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
remove_diamond_seed_genes <- function(diamond_module){
  diamond_module$module_genes <- setdiff(diamond_module$module_genes, diamond_module$seed_genes)
  diamond_module$settings$include_seed <- F
  return(diamond_module)
}
#' Add seed genes from a DIAMOnD MODifieR_module
#' @inheritParams remove_diamond_seed_genes
#' @details 
#' Adds seed genes from a DIAMOnD module
#' @return
#' An object of class "MODifieR_module" with subclass "DIAMOnD"
#' @seealso 
#' \code{\link{diamond}}
#' @export
#' @export
add_diamond_seed_genes <- function(diamond_module){
  diamond_module$module_genes <- unique(c(diamond_module$module_genes, diamond_module$seed_genes))
  diamond_module$settings$include_seed <- T
  return(diamond_module)
}