#' MCODE
#' 
#' A clique based Algorithm to identify disease module from differentially expressed genes 
#' 
#' @param deg_genes A nx2 dataframe sorted ranked list of genes with their respective -log10 p-values 
#' @param type keytype of the gene that is provide (default = "Gene Symbol")
#' @param hierarchy the hierarchical level of the PPI network to be used (default = 1)
#' @param vwp Vertex weight percentage. Default value is 0.5.
#' @param haircut	Boolean value, whether to remove singly-connected nodes from clusters (TRUE) or not (FALSE).
#' @param fluff	Boolean value, whether to spand cluster cores by one neighbour shell outwards (TRUE) or not (FALSE).
#' @param fdt	Cluster density cutoff. Default value is 0.8.
#' @param loops	Boolean value, whether to include self-loops (TRUE) or not (FALSE).
#' @export
mod_mcode <- function(MODifieR_input, ppi_network, type = "Gene symbol", hierarchy = 1, vwp =0.5, haircut = F, fluff = F,
                  fdt = 0.8, loops = T, diffgen_cutoff = 0.05, module_cutoff = 3.5, dataset_name = NULL){
  default_args <- formals()
  user_args <- as.list(match.call(expand.dots = T)[-1])
  settings <- c(user_args, default_args[!names(default_args) %in% names(user_args)])

  if (!is.null(dataset_name)){
    settings$MODifieR_input <- dataset_name
  }

  deg_genes <- MODifieR_input$diff_genes
  
  deg_genes$pvalue <- -log10(deg_genes$pvalue)
  
  deg_genes <- deg_genes[order(deg_genes$pvalue, decreasing = T),]
  
  colnames(deg_genes) <- c("hgnc_symbol" , "p_val")
  deg_genes <- deg_genes[deg_genes$p_val > (-log10(diffgen_cutoff)), ]
 
  colnames(ppi_network) <- c("Interactor.1.Gene.symbol", "Interactor.2.Gene.symbol")
  
  network <-construction_mod(input = as.data.frame(deg_genes$hgnc_symbol),
                         db = ppi_network,
                         species = "human",
                         ID.type = type,
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
