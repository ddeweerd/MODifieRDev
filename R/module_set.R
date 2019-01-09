#'Consensus module
#'
#'@param min_frequency Minimal number of MODifieR modules that a gene should be present in in order
#'to include it in the final module
#'@param module_list A list of MODifieR modules
#'@details
#'Get a consensus module that is composed of genes present in at least \code{min_frequency} genes.
#'If the input \code{module_list} is unnamed, the subclasses of the MODifieR objects will be used.
#'If there is more than 1 of a given subclass present in \code{module_list}, a number will be appended.
#'
#'@return
#' Returns an object of class "MODifieR_module" with subclass "Module_set". 
#' This object is a named list containing the following components:
#' \item{module_genes}{A character vector containing the genes in the final module}
#' \item{module_gene_list}{A named list containing the module genes from the original modules}
#' \item{gene_frequency}{Table containing all the genes present in the modules and their frequency}
#' \item{method_by_gene}{A named list where the elements are the modules the genes have been found in
#' and the names are the gene names}
#' \item{gene_by_method}{A table containing the gene frequencies by combination of methods}
#' \item{settings}{A named list containing the parameters used in generating the object}
#'@export
create_module_set <- function(min_frequency, module_list){
  
  # Retrieve settings
  settings <- do.call(what = "settings_function", as.list(stackoverflow::match.call.defaults()[-1]))
  
  module_gene_list <- sapply(X = module_list, FUN = function(x)x$module_genes)
  
  if (is.null(names(module_list))){
    names(module_gene_list) <- sapply(X = module_list, FUN = function(x){class(x)[2]})
  }else{
    names(module_list)[sapply(names(module_list), function(x)x=="")] <- 
      sapply(X = module_list[sapply(names(module_list), function(x)x=="")] , FUN = function(x){class(x)[2]})
    names(module_gene_list) <- names(module_list)
  }
  
  tabled_names <- table(names(module_gene_list))
  
  duplicated_names <- names(tabled_names[tabled_names > 1])
  
  if (sum(tabled_names[tabled_names > 1]) != 0){
    for (i in 1:length(duplicated_names)){
      message(duplicated_names[i], " is non-unique, appending number to ", duplicated_names[i])
      names(module_gene_list)[grep(pattern = duplicated_names[i], x = names(module_gene_list))] <- 
        paste0(names(module_gene_list)[grep(pattern = duplicated_names[i], x = names(module_gene_list))], 1:tabled_names[tabled_names > 1][i])
    }
  }
  
  gene_frequency <- table(unlist(module_gene_list))
  
  if (max(gene_frequency) < min_frequency){
    stop("Minimum frequency higher than maximal overlap in modules. Maximal overlap is ", max(gene_frequency))
  }
  module_genes <- names(gene_frequency[gene_frequency >= min_frequency])
  
  gene_table <- table_gene_by_method(genes = module_genes, module_gene_list = module_gene_list)
  
  method_by_gene <- gene_table[[1]]
  names(method_by_gene) <- module_genes
  
  gene_by_method = gene_table[[2]]
  
  new_module_set <- construct_module_set(module_genes = module_genes, 
                                         module_gene_list = module_gene_list,
                                         gene_frequency = gene_frequency,
                                         method_by_gene = method_by_gene,
                                         gene_by_method = gene_by_method,
                                         settings = settings)
  return(new_module_set)
}

construct_module_set <- function(module_genes, module_gene_list, gene_frequency, method_by_gene, gene_by_method, settings){
  new_module_set <- list("module_genes" = module_genes,
                         "module_gene_list"= module_gene_list,
                         "gene_frequency" = gene_frequency,
                         "method_by_gene" = method_by_gene,
                         "gene_by_method" = gene_by_method,
                         "settings" = settings)
  
  class(new_module_set) <- c("MODifieR_module", "Module_set")
  
  return(new_module_set)
}

table_gene_by_method <- function(genes, module_gene_list){
  gene_frequencies <- lapply(genes, FUN = get_gene_presence, module_gene_list = module_gene_list)
  gene_table <- get_gene_table(gene_frequencies = gene_frequencies)
  return(list("gene_frequencies" = gene_frequencies, "module_frequency" = gene_table))
}

get_gene_presence <- function(gene, module_gene_list){
  names(module_gene_list)[sapply(1:length(module_gene_list), function(x)gene %in% module_gene_list[[x]])]
}
#' Get maximal number of MODifieR modules that share at least one gene
#' 
#' @inheritParams create_module_set
#' 
#' @return
#' The function returns an integer that gives the maximal number of MODifieR objects in 
#' \code{module_list} that share at least one gene.
#' @export
get_max_frequency <- function(module_list){
  max(table(unlist(sapply(module_list, function(x)x$module_genes))))
}

get_gene_table <- function(gene_frequencies){
  max_combinations <- max(lengths(gene_frequencies))
  names_combinations <- unique(unlist(gene_frequencies))
  
  results <- lapply(X = 2:max_combinations, FUN = get_gene_frequencies, 
                    names_combinations = names_combinations, 
                    gene_frequencies = gene_frequencies)
  
  results <- unlist(results)
  
  results <- results[results  > 0 ]
  
  return (results)
}
get_gene_frequencies <- function(m, names_combinations, gene_frequencies){
  combination_table <- get_combinations(m = m, x = names_combinations)
  frequencies <- apply(X = combination_table, MARGIN = 2, get_intersections, gene_frequencies = gene_frequencies)
  names(frequencies) <- apply(combination_table, MARGIN = 2, FUN = paste, collapse = " ")
  
  return (frequencies)
}

get_combinations <- function(m, x){
  combn(x = x, m = m)
}

get_intersections <- function(combo, gene_frequencies){
  sum(sapply(X = gene_frequencies, FUN = function(x){sum(combo %in% x) == length(combo)}))
}
