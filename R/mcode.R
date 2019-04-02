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
#' 
#' @export
mod_mcode <- function(MODifieR_input, ppi_network, hierarchy = 1, vwp =0.5, haircut = F, fluff = F,
                  fdt = 0.8, loops = T, deg_cutoff = 0.05, module_cutoff = 3.5, dataset_name = NULL){
 
  # Retrieve settings
  evaluated_args <- c(as.list(environment()))
  settings <- as.list(stackoverflow::match.call.defaults()[-1])
  replace_args <- names(settings)[!names(settings) %in% unevaluated_args]
  for (argument in replace_args) {
    settings[[which(names(settings) == argument)]] <- evaluated_args[[which(names(evaluated_args) == 
                                                                              argument)]]
  }
  
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
                         hierarchy = 1)

  result <- mcode(network, vwp = vwp, haircut = haircut,
                          fluff = fluff, fdt = fdt)

  result_genes <- result$COMPLEX[result$score > module_cutoff]
  result_scores <- result$score[result$score > module_cutoff]
  
  result_df <- c()
  for (k in 1:length(result_genes)){
    result_df[[k]] <- ppi_network[result_genes[[k]] ,]
    
  }
  

  result_modules <- list()

  for (i in 1:length(result_genes)){
    result_modules[[i]] <- mcode_module(module_genes =  unique(c(as.character(result_df[[i]]$Interactor.1.Gene.symbol) , as.character(result_df[[i]]$Interactor.2.Gene.symbol))),
                                        module_score = result_scores[i],
                                        settings = settings)
  }

  return(result_modules)
}

mcode_module <- function(module_genes, module_score, settings){
  new_mcode_module <- list("module_genes" = module_genes,
                           "module_score" = module_score,
                           "settings " = settings)

  class(new_mcode_module) <- c("MODifieR_module", "Mcode")
  return(new_mcode_module)
}
