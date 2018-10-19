#' @import org.Hs.eg.db
#' @import igraph
#' @import foreach
#' @import doParallel
#' @import WGCNA
#' @importFrom Rcpp sourceCpp
#' @importFrom flashClust flashClust
#' @useDynLib MODifieRDev
#' @importFrom dynamicTreeCut printFlush

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
  conv <- toTable(org.Hs.eg.db::org.Hs.egSYMBOL)
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
  conv <- toTable(org.Hs.eg.db::org.Hs.egSYMBOL)
  MODifieR_module$module_genes <- conv$symbol[conv$gene_id %in% MODifieR_module$module_genes]
  return(MODifieR_module)
}
dataframe_to_vector <- function(dataframe){
  setNames(dataframe[,2],dataframe[,1])
}

summary.MODifieR_module <- function(MODifieR_module){
 length(MODifieR_module$module_genes)
}

settings_function <- function(...) {
  func_args <- sapply(match.call(expand.dots=TRUE)[-1], deparse)
  func_args <- gsub(pattern = "\"", replacement = "", x = func_args)
  return(func_args)
}
#' @export
get_string_DB_ppi <- function(version = "10", score_threshold = 900, simplify_graph = T){
  string_db <- STRINGdb$new( version=as.character(version), species=9606,
                             score_threshold=score_threshold, input_directory="" )
  ppi_network <- string_db$get_graph()

  if (simplify_graph == T){
    ppi_network <- igraph::simplify(ppi_network)
  }

  return (ppi_network)
}

#' @export
symbol_to_prot <- function(dataframe, remove_NA = T){
  prot_ids_list <-  mapIds(org.Hs.eg.db,
                      keys=dataframe[,1],
                      column="ENSEMBLPROT",
                      keytype="SYMBOL",
                      multiVals="list")
  prot_ids_unlisted <- lapply(X = 1:length(prot_ids_list), FUN = function(x){
    cbind(names(prot_ids_list[x]), unlist(prot_ids_list[x], use.names = F))
  })

  prot_ids_final <- base::Reduce(f = rbind, x = prot_ids_unlisted)
  if (remove_NA == T){
    prot_ids_final <- na.omit(prot_ids_final)
  }

  colnames(prot_ids_final) <- c("hgnc_symbol", "ensembl_peptide_id")

  return(prot_ids_final)
}
#' @export
prot_to_entrez <- function(dataframe){
  entrez_ids <-  na.omit(t(as.data.frame(mapIds(org.Hs.eg.db,
                        keys= dataframe[,1],
                        column="ENTREZID",
                        keytype="ENSEMBLPROT",
                        multiVals="list"))))
  return(entrez_ids[,1])
}

#' @export
extract_module_data <- function(MODifieR_module, data_field){
  get(data_field, MODifieR_module)
}
#' @export
extract_module_class <- function(MODifieR_module){
  class(MODifieR_module)[2]
}
#' @export
#Unsure if this should be included in the package
print_super_module <- function(super_module, output_file){
  sapply(X = 1:length(super_module), FUN =  function(i){
  invisible(write.table(x = t(c(names(super_module[i]), "", super_module[[i]])),
                  file = output_file, append = T, row.names = F,
                  col.names = F, sep = "\t", quote = F))})
}
