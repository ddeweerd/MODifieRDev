#' @import org.Hs.eg.db
#' @import igraph
#' @import foreach
#' @import doParallel
#' @import parallel
#' @import WGCNA
#' @import reticulate
#' @import MODA
#' @import AnnotationDbi
#' @import RSQLite
#' @import edgeR
#' @importFrom plyr ddply summarise
#' @importFrom stats as.dist fisher.test hclust median na.omit p.adjust phyper pnorm
#' @importFrom stats prcomp pt qt quantile runif var
#' @importFrom Rcpp sourceCpp
#' @importFrom flashClust flashClust
#' @importFrom dynamicTreeCut printFlush
#' @importFrom utils combn
#' @importFrom graphics plot
#' @importFrom DESeq2 varianceStabilizingTransformation
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom stackoverflow match.call.defaults
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


get_diffgene_cutoff <- function(gene_fraction, MODifieR_module){
  sort(MODifieR_module$diff_genes$pvalue)[round(length(MODifieR_module$diff_genes$gene) * gene_fraction)]
}



settings_function <- function(...){
  evaluated_args <- c(as.list(environment()), list(...))
  settings <- as.list(stackoverflow::match.call.defaults()[-1])
  replace_args <- names(settings)[!names(settings) %in% unevaluated_args]
  for (argument in replace_args) {
    settings[[which(names(settings) == argument)]] <- evaluated_args[[which(names(evaluated_args) == 
                                                                              argument)]]
  }
  return(settings)
}
#Generic subclass extractor
extract_module_class <- function(MODifieR_module){
  class(MODifieR_module)[2]
}
#'@export
plot.MODifieR_module <-function(MODifieR_module, ppi_network){
  ppi_graphed <- module_to_igraph(MODifieR_module = MODifieR_module, ppi_network = ppi_network)
  deg <- degree(ppi_graphed, mode="all")
  V(ppi_graphed)$size <- deg /10
  graphics::plot(ppi_graphed, edge.arrow.size=.1, vertex.label.cex=.45)
}

module_to_igraph <- function(MODifieR_module, ppi_network){
  genes <- MODifieR_module$module_genes[MODifieR_module$module_genes %in% unique(unname(unlist(ppi_network[,1:2])))]
  ppi_graphed <- igraph::graph.data.frame(ppi_network)
  ppi_graphed <- igraph::simplify(ppi_graphed, edge.attr.comb=list(Weight="sum","ignore"))
  ppi_graphed <- igraph::induced.subgraph(graph = ppi_graphed, vids = genes)
}

summary.MODifieR_input <- function(MODifieR_input){
  MODifieR_input$diff_genes
  
  significance_levels <- c(0.05, 0.01, 0.001, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10)
  
  n_of_genes <- sapply(X = significance_levels, FUN = function(x)sum(MODifieR_input$diff_genes$pvalue < x))
  
  names(n_of_genes) <- significance_levels
  
  return(list(significant_genes = n_of_genes))
  
}

