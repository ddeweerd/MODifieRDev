#' MODA
#'
#' An implementation of MODA co-expression based algorithm.
#' @inheritParams clique_sum
#' @inheritParams MODA::CompareAllNets
#' @inheritParams MODA::WeightedModulePartitionHierarchical
#' @param group_of_interest Numerical value denoting which group contains the condition of interest (1 or 2) 
#' 
#' @details 
#' This implementation follows a workflow as described in the MODA vignette. First,
#' two separate networks are constructed, a background network containing expression
#' data from all samples and a condition specific network consisting of all samples minus
#' the condition specific samples.
#' Then, hierarchical clustering is performed and cutting height estimated from either 
#' maximal average \code{density} or \code{modularity}
#' 
#' \emph{Condition} specific co-expression modules are then extracted 
#' using the Jaccard index and \code{specificTheta}. 
#' 
#' The final module will consist of the co-expression module that has the minimal 
#' Jaccard index complemented by co-expression modules that have a Jaccard index 
#' below this minimal + \code{specificTheta}
#' 
#' After analysis, the \code{specificTheta}
#' and thereby the disease module can be adjusted using 
#' \code{\link{moda_change_specific_threshold}}
#' 
#' @return moda returns an object of class "MODifieR_module" with subclass "MODA". 
#' This object is a named list containing the following components:
#' \item{module_genes}{A character vector containing the genes in the final module}
#' \item{group1_modules}{A list containing all co-expression modules in the background network}
#' \item{group2_modules}{A list containing all co-expression modules in the condition specific network}
#' \item{jaccard_table}{A matrix with all Jaccard indexes for all co-expression modules}
#' \item{settings}{A named list containing the parameters used in generating the object}
#' 
#' @seealso 
#' \url{https://bioconductor.org/packages/release/bioc/vignettes/MODA/inst/doc/MODA.html}
#' @references 
#' \cite{Li D, Brown JB, Orsini L, Pan Z, Hu G, He S (2016). MODA: MODA: 
#' MOdule Differential Analysis for weighted gene co-expression network. R package version 1.6.0}
#' @export
moda <- function(MODifieR_input,
                   cutmethod = "Density",
                   specificTheta = 0.1, conservedTheta = 0.1,
                   dataset_name = NULL, group_of_interest){
  # Retrieve settings
  settings <- do.call(what = "settings_function", as.list(stackoverflow::match.call.defaults()[-1]))
  
  if (!is.null(dataset_name)){
    settings$MODifieR_input <- dataset_name
  }
  ##Get relevant input data from input object 
  #Background network, meaning data from all conditions
  datExpr1 <- t(MODifieR_input$annotated_exprs_matrix)
  #Condition specific network, meaning all samples MINUS the samples belonging to the condition
  datExpr2 <- t(MODifieR_input$annotated_exprs_matrix[,-MODifieR_input$group_indici[[group_of_interest]]])
  
  indicator1 = names(MODifieR_input$group_indici)[1]  # indicator for data profile 1
  indicator2 = names(MODifieR_input$group_indici)[2]  # indicator for data profile 2
  
  modules1 <- WMPH(datExpr = datExpr1, indicatename = indicator1, cutmethod = cutmethod)
  
  modules2 <- WMPH(datExpr = datExpr2, indicatename = indicator2, cutmethod = cutmethod)
  
  
  jaccard_table <- compare_modules(module_list1 = modules1, module_list2 = modules2)
  
  module_genes <- unique(unlist(modules1[moda_extract_modules_index_specific(
    jaccard_table = jaccard_table, specificTheta = specificTheta)]))
  
  
  
  new_moda_module <- list("module_genes" =  module_genes,
                          "group1_modules" = modules1,
                          "group2_modules" = modules2,
                          "jaccard_table" = jaccard_table,
                          "settings" = settings)
  
  class(new_moda_module) <- c("MODifieR_module", "MODA")
  
  return (new_moda_module)
}
#'
#'Change 
#'@inheritParams MODA::CompareAllNets
#'@param moda_module  A \code{MODifieR_input} object created by \code{\link{moda}}
#'
#'@export
moda_change_specific_threshold <- function(moda_module, specificTheta){
  
  module_genes <- unique(unlist(moda_module$group1_modules[moda_extract_modules_index_specific(
    jaccard_table = moda_module$jaccard_table, specificTheta = specificTheta)]))
  
  moda_module$module_genes <- module_genes
  moda_module$settings$specificTheta <- specificTheta
  
  return(moda_module)
}
