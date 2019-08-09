#' MCODE
#' 
#' A clique based algorithm to identify disease modules from differentially expressed genes 
#' originally by Bader et al
#' 
#' @details 
#' 
#' Much of the code an documentation has been taken from the now defunct package "ProNet"
#' 
#' @inheritParams clique_sum_exact
#' @inheritParams build_clique_db
#' @param hierarchy This parameter indicates how many hierarchy are included in the network, currently it can be \code{0}, \code{1} or \code{2}. Default value is \code{1}.
#' @param vwp Vertex weight percentage. Default value is 0.5.
#' @param haircut	Boolean value, whether to remove singly-connected nodes from clusters (TRUE) or not (FALSE).
#' @param fluff	Boolean value, whether to spand cluster cores by one neighbour shell outwards (TRUE) or not (FALSE).
#' @param fdt	Cluster density cutoff. Default value is 0.8.
#' @param loops	Boolean value, whether to include self-loops (TRUE) or not (FALSE).
#' @param module_cutoff Minimal score for a module to be returned
#' 
#' @return 
#' mcode returns a list of objects of class "MODifieR_module" with subclass "Mcode". 
#' The objects are named lists containing the following components:
#' \item{module_genes}{A character vector containing the genes in the final module}
#' \item{module_scores}{A numeric value that denotes the score of the module. Higher is better}
#' \item{settings}{A named list containing the parameters used in generating the object}
#' 
#' @seealso 
#' 
#' \url{https://github.com/cran/ProNet}
#' 
#' @references \cite{Bader, G. D., & Hogue, C. W. (2003). An automated method for finding molecular complexes in 
#' large protein interaction networks. BMC Bioinformatics, 4(1), 2. \url{https://doi.org/10.1186/1471-2105-4-2}}
#' @author DIrk de Weerd
#' 
#' @export
mod_mcode <- function(MODifieR_input, ppi_network, hierarchy = 1, vwp = 0.5, haircut = F, fluff = F,
                      fdt = 0.8, loops = T, deg_cutoff = 0.05, module_cutoff = 3.5, dataset_name = NULL){
  
  # Retrieve settings
  evaluated_args <- c(as.list(environment()))
  settings <- as.list(stackoverflow::match.call.defaults()[-1])
  replace_args <- names(settings)[!names(settings) %in% unevaluated_args]
  for (argument in replace_args) {
    settings[[which(names(settings) == argument)]] <- evaluated_args[[which(names(evaluated_args) == 
                                                                              argument)]]
  }
  
  #Validate the input parameters
  check_diff_genes(MODifieR_input, deg_cutoff = deg_cutoff)
  ppi_network <- validate_ppi(ppi_network)
  validate_inputs(settings)
  
  if (!is.null(dataset_name)){
    settings$MODifieR_input <- dataset_name
  }
  
  deg_genes <- MODifieR_input$diff_genes
  
  deg_genes$pvalue <- -log10(deg_genes$pvalue)
  
  deg_genes <- deg_genes[order(deg_genes$pvalue, decreasing = T),]
  
  colnames(deg_genes) <- c("hgnc_symbol" , "p_val")
  deg_genes <- deg_genes[deg_genes$p_val > (-log10(deg_cutoff)), ]
  
  colnames(ppi_network) <- c("Interactor.1.Gene.symbol", "Interactor.2.Gene.symbol")
  
  network <-construction_mod(input = as.data.frame(deg_genes$hgnc_symbol),
                             db = ppi_network,
                             species = "human",
                             hierarchy = hierarchy)
  
  result <- mcode(network, vwp = vwp, haircut = haircut,
                  fluff = fluff, fdt = fdt)
  
  result_genes <- result$COMPLEX[result$score > module_cutoff]
  result_scores <- result$score[result$score > module_cutoff]

  result_df <- lapply(X = result$COMPLEX, FUN = index_to_gene, ppi_network = ppi_network)
  result_df <- lapply(result_df, get_module_genes)
  
  module_genes <- get_module_genes(result_df[result$score > module_cutoff])
   
  result_modules <- construct_mcode_module(module_genes =  module_genes,
                                      modules = result_df,
                                      module_scores = result$score,
                                      input_class = class(MODifieR_input)[3],
                                      settings = settings)
  
  
  return(result_modules)
}
#@author Dirk de Weerd
construct_mcode_module <- function(module_genes, modules, module_scores, input_class, settings){
  new_mcode_module <- list("module_genes" = module_genes,
                           "modules" = modules,
                           "module_scores" = module_scores,
                           "settings" = settings)
  
  class(new_mcode_module) <- c("MODifieR_module", "Mcode", input_class)
  return(new_mcode_module)
}

index_to_gene <- function(indici, ppi_network){
  ppi_network[indici, ]
}

get_module_genes <- function(gene_set){
 unique(unname(unlist(gene_set)))
}
#' @export
mcode_adjust_score <- function(module_cutoff, mcode_module){
  mcode_module$settings$module_cutoff <- module_cutoff
  
  module_genes <-  get_module_genes(mcode_module$modules[mcode_module$module_score > module_cutoff])
  
  new_mcode_module <- construct_mcode_module(module_genes = module_genes, 
                                   modules = mcode_module$modules,
                                   module_scores = mcode_module$module_scores,
                                   input_class = class(mcode_module)[3],
                                   settings = mcode_module$settings)
  
  return(new_mcode_module)
}
#' @export
mcode_split_modules <- function(module_cutoff, mcode_module){
  new_mcode_modules <- list()
  for(i in 1:length(mcode_module$modules[mcode_module$module_scores > module_cutoff])){
    settings <- mcode_module$settings
    settings$module_cutoff <- mcode_module$module_scores[i]
    new_mcode_modules[[i]] <- construct_mcode_module(module_genes = mcode_module$modules[[i]], 
                                           modules = , mcode_module$modules, 
                                           module_scores <- mcode_module$module_scores,
                                           input_class = class(mcode_module)[3],
                                           settings = settings)
  }
  return(new_mcode_modules)
}