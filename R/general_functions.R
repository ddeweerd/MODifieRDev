#' @import org.Hs.eg.db
#' @import igraph
#' @import foreach
#' @import doParallel
#' @import parallel
#' @import WGCNA
#' @import reticulate
#' @import MODA
#' @import AnnotationDbi
#' @importFrom plyr ddply summarise
#' @importFrom stats as.dist fisher.test hclust median na.omit p.adjust phyper pnorm
#' @importFrom stats prcomp pt qt quantile runif
#' @importFrom Rcpp sourceCpp
#' @importFrom flashClust flashClust
#' @importFrom dynamicTreeCut printFlush
#' @importFrom utils combn
#' @useDynLib MODifieRDev

#'@title Convert the module genes in a MODifieR_input object from official gene symbols to ENTREZ gene IDs
#' @param MODifieR_module An object of class MODifieR_module
#' @details 
#' The function uses the \code{org.Hs.egSYMBOL} function from the package \code{org.Hs.eg.db}
#' to convert ENTREZ IDs to official gene symbols
#' @seealso 
#' \code{\link[org.Hs.eg.db]{org.Hs.egSYMBOL}}
#' 
#'  \code{\link{entrez_to_symbol}}
#' @export
symbol_to_entrez <- function(MODifieR_module){
  conv <- AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egSYMBOL)
  MODifieR_module$module_genes <- conv$gene_id[conv$symbol %in% MODifieR_module$module_genes]
  return(MODifieR_module)
}

#' Convert the module genes in a MODifieR_input object from ENTREZ gene IDs to official gene symbols
#' @param MODifieR_module An object of class MODifieR_module
#' @details 
#' The function uses the \code{org.Hs.egSYMBOL} function from the package \code{org.Hs.eg.db}
#' to convert official gene symbols to ENTREZ IDs
#' @seealso 
#' \code{\link[org.Hs.eg.db]{org.Hs.egSYMBOL}}
#' 
#' \code{\link{symbol_to_entrez}}
#' @export
entrez_to_symbol <- function(MODifieR_module){
  conv <- AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egSYMBOL)
  MODifieR_module$module_genes <- conv$symbol[conv$gene_id %in% MODifieR_module$module_genes]
  return(MODifieR_module)
}
dataframe_to_vector <- function(dataframe){
  stats::setNames(dataframe[,2],dataframe[,1])
}

summary.MODifieR_module <- function(MODifieR_module){
 length(MODifieR_module$module_genes)
}

settings_function <- function(...) {
  func_args <- sapply(match.call(expand.dots=TRUE)[-1], deparse)
  func_args <- gsub(pattern = "\"", replacement = "", x = func_args)
  return(func_args)
}
#Generic subclass extractor
extract_module_class <- function(MODifieR_module){
  class(MODifieR_module)[2]
}
