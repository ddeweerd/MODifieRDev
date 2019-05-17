#' Clique_correlation
#'
#' A clique based method to find a disease module from correlated gene expression
#'
#' @inheritParams clique_sum_exact
#' @inheritParams build_clique_db
#' @param iteration Number of iterations to be performed
#' @param fraction_of_interactions Fraction of interactions from the original network that will be used
#' in each iteration
#' @param frequency_cutoff Fraction of the number of times a gene should be present in it 
#' iterations. Default is 0.5, meaning 50 procent of all iterations 
#' @param clique_significance Cutoff for Fisher exact test for cliques
#' @param multiple_cores parallelize iterations using number of cores on system -1?
#' @param tempfolder Folder where temporary files are stored
#' @details 
#' 
#' The correlation clique is a clique-based algorithm using consensus clustering.
#' The algorithm starts with calculating a correlation score between each interaction in the PPi network.
#' The correlation score is obtained by subtracting the Pearson correlation p-value:
#' \deqn{correlation score = 1 - Pearson p-value}
#' Subsequently, the correlation score is multiplied by 
#' the correlation confidence and scaled with the scale factor to get the edge score:
#' \deqn{edge score = \sqrt(correlation score * confidence score  ) * scale_factor}
#' When the edge scores are calculated the iterative part of the algorithm commences:
#' All edge scores are compared to random variables from the uniform distribution between (0,1)
#' only interactions where the edge score is higher than the random variable are used to construct
#' a new PPi network. Then, maximal cliques are inferred from this new network.
#' The cliques are tested for significant enrichment of DEGs by Fisher's exact test and
#' the union of significant cliques is the disease module for this iteration
#' The final disease module will consist of genes that have been present in at least \code{frequency_cutoff} iterations
#' 
#' @return 
#' correlation_clique returns an object of class "MODifieR_module" with subclass "Correlation_clique". 
#' This object is a named list containing the following components:
#' \item{module_genes}{A character vector containing the genes in the final module}
#' \item{frequency_table}{A table containing the fraction of times the genes were present in an iteration module}
#' \item{settings}{A named list containing the parameters used in generating the object}
#' 
#' @author Dirk de Weerd
#' 
#' @export
correlation_clique <- function(MODifieR_input, ppi_network,
                               frequency_cutoff = .5, 
                               fraction_of_interactions = 0.4,
                               iteration = 50, clique_significance = 0.01, 
                               deg_cutoff = 0.05, multiple_cores = F, n_cores = 3, 
                               dataset_name = NULL){
  
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
  MODifieR_input$diff_genes <- MODifieR_input$diff_genes[MODifieR_input$diff_genes$gene %in% rownames(MODifieR_input$annotated_exprs_matrix), ]
  #Making sure the PPi network has 3 columns
  springConnection <- ppi_network[,1:3]
  #Get DEGs
  pValueMatrix <- MODifieR_input$diff_genes
  pValueMatrix <- stats::na.omit(pValueMatrix)
  #Remove genes from PPi that are not in input
  overlap <- springConnection[,1] %in% pValueMatrix[,1] & springConnection[,1] %in% rownames(MODifieR_input$annotated_exprs_matrix)
  springConnection <- springConnection[overlap,]
  overlap <- springConnection[,2] %in% pValueMatrix[,1] & springConnection[,2] %in% rownames(MODifieR_input$annotated_exprs_matrix)
  springConnection <- springConnection[overlap,]
  
  springConnection[,3] <- springConnection[,3] / 1000
  #Calculates (1- p value) for all relevant connections in the PPi network, and is stored in column 1.
  #Column 2 will contain the PPi score from the PPi network
  corrPvalues <- cbind(apply(X = springConnection, MARGIN = 1, FUN = calculate_correlation, 
                             expression_matrix = MODifieR_input$annotated_exprs_matrix),
                       springConnection[,3])
  #Calculates square root of p-value * PPi score. This is then scaled by probabilityScaleFactor. 
  pval_score <- apply(X = corrPvalues, MARGIN = 1, FUN = function(x){sqrt(x[1] * x[2])})
  
  pval_score <- scale_observations(sf = fraction_of_interactions, scores = pval_score)
  
  genes <- dataframe_to_vector(as.data.frame(pValueMatrix))
  
  genes <- genes[names(genes) %in% unique(unlist(springConnection[ ,1:2]))]
  
  module_list <- list()
  if (multiple_cores){
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    parallel::clusterCall(cl, function(x) .libPaths(x), .libPaths())
    #module_list <- foreach(i=1:iteration, .combine = list, .export = ls(.GlobalEnv)) %do% {
    #Compare computed scores to random scores
    
    module_list <- parLapply(cl = cl, X = 1:iteration, fun = perform_iterations, pval_score, springConnection, genes, clique_significance, deg_cutoff)
    parallel::stopCluster(cl) # stop the cluster
    
    #}
    
  }else{
    #For every iteration:
    for (i in 1:iteration){
      #Compare computed scores to random scores 
      
      
      module_list[[i]] <- perform_iterations(iteration_n = i, pval_score = pval_score, springConnection = springConnection, genes = genes, clique_significance = clique_significance, deg_cutoff = deg_cutoff)
    }
  }
  #Get the frequency of each gene that has been present in a module, expressed in fraction
  tabled_frequencies <- table(unlist(module_list)) / iteration
  #The resulting module is going to consist of all genes that have are seen more frequently in the iterations
  #than the frequency cutoff (default 0.5, so present in 50% of all iterations) 
  module_genes <- names(tabled_frequencies[tabled_frequencies >= frequency_cutoff])
  
  new_correlation_clique_module <- construct_correlation_module(module_genes = module_genes,
                                                                frequency_table = tabled_frequencies,
                                                                settings = settings)
  
  return( new_correlation_clique_module)
}


perform_iterations <- function(iteration_n, pval_score, springConnection, genes, clique_significance, deg_cutoff){
  
  to_pass <- pval_score > runif(n = length(pval_score))
  
  clique_file_list <- create_clique_file(adjacency_matrix = springConnection[to_pass, ], genes = genes) 
  
  module_list <- process_clique(clique_file = clique_file_list$tempfile, genes = clique_file_list$tr_genes, 
                                graphed_ppi = clique_file_list$graphed_ppi, clique_significance = clique_significance,
                                deg_cutoff = deg_cutoff, name_table = clique_file_list$name_table)
}
#Pval Helper function
fisher_pval <- function(clique_row){
  stats::fisher.test(x = matrix(data = clique_row, nrow = 2, byrow = T), alternative = "g")$p.value
}
#This function takes as the first argument a row from a PPi network, retrieves the corresponding row and column
#by row and column name and computes (1- correlation P value). 
calculate_correlation <- function(row, expression_matrix){
  row <- as.character(row)
  
  x_y <- c(which(rownames(expression_matrix) == row[1]), which(rownames(expression_matrix) == row[2]))
  
  1 - round(stats::cor.test(x = expression_matrix[x_y[1],], y = expression_matrix[x_y[2],])$p.value, 3)
}

#' correlation_adjust_cutoff
#' 
#' @inheritParams correlation_clique
#' @param correlation_module Module object that has been produced by \code{correlation_clique}  function
#' @details 
#' 
#' This function allows to adjust the frequency cutoff for a \code{correlation_clique} module object
#' @return 
#' 
#' 
#'  \code{correlation_clique} module object 
#' @seealso 
#' 
#' \code{\link{correlation_clique}}
#' 
#' @author Dirk de Weerd
#' 
#' @export
correlation_adjust_cutoff <- function(frequency_cutoff, correlation_module){
  
  frequency_table <- correlation_module$frequency_table
  
  module_genes <- names(frequency_table[frequency_table >= frequency_cutoff])
  
  correlation_module$settings$frequency_cutoff <- frequency_cutoff
  
  new_correlation_clique_module <- construct_correlation_module(module_genes = module_genes,
                                                                frequency_table = frequency_table,
                                                                settings = correlation_module$settings)
  return( new_correlation_clique_module)
  
}

#' correlation_set_module_size
#' 
#' Returns a correlation_clique module closest to \code{size}
#' 
#' @inheritParams correlation_adjust_cutoff
#' @param size Module object that has been produced by \code{correlation_clique}
#'  function
#' @details 
#' The function will find the the frequency cutoff for that will result
#' in a \code{correlation_clique} module object closest to \code{size}
#' @return 
#'  \code{correlation_clique} module object 
#' @seealso 
#' 
#' \code{\link{correlation_clique}}
#' 
#' @author Dirk de Weerd
#' 
#' @export
correlation_set_module_size <- function(size, correlation_module){
  frequency_cutoff <- as.numeric(names(which.min(abs(size - table(correlation_module$frequency_table)))))
  
  new_correlation_clique_module <- correlation_adjust_cutoff(frequency_cutoff = frequency_cutoff,
                                                             correlation_module = correlation_module)
  
  return (new_correlation_clique_module)
}

construct_correlation_module <- function(module_genes, frequency_table, settings){
  
  new_correlation_clique_module <- list("module_genes" = module_genes,
                                        "frequency_table" = frequency_table,
                                        "settings" = settings)
  
  class( new_correlation_clique_module) <- c("MODifieR_module", "Correlation_clique")
  
  return( new_correlation_clique_module)
  
}

create_clique_file <- function(adjacency_matrix, genes){
  clique_file <- tempfile()
  
  g <- igraph::simplify(igraph::graph.data.frame(adjacency_matrix, directed = F))
  
  genes <- genes[names(genes) %in% V(g)$name]
  name_table <- cbind(V(g)$name, as.character(as.numeric(V(g)-1)))
  
  names(genes) <- sapply(X = names(genes), FUN = entrez_to_index, name_table = name_table)
  
  igraph::max_cliques(graph = g, file = clique_file)

  return(list(tempfile = clique_file, graphed_ppi = g, name_table = name_table, tr_genes = genes))
}


process_clique <- function(clique_file, genes, graphed_ppi, clique_significance, deg_cutoff, name_table) {
 
  deg_genes <- genes[genes < deg_cutoff]
  
  non_deg_genes <- genes[genes >= deg_cutoff]
  cutoff_per_size <- vector(mode = "numeric")
  con <- file(clique_file, "r")
  significant_cliques <- NULL
  deg_names <- names(deg_genes)
  while (TRUE) {
    line <- readLines(con, n = 100000)
    if ( length(line) == 0 ) {
      break
    }
    clique_genes <- sapply(X = line, FUN = strsplit, split = " ")
    
    system.time(size_and_deg <- (unname(sapply(X = clique_genes, FUN = function(x) c(length(x), sum(x %in% deg_names))))))
    
    clique_sizes <- unique(size_and_deg[1,])
    
    if (sum(is.na(cutoff_per_size[clique_sizes]) > 0)){
      indici <- clique_sizes[which(is.na(cutoff_per_size[clique_sizes]))]
      cutoff_per_size[indici] <- sapply(X = indici, FUN = get_clique_cutoffs, 
                                        deg_names = deg_names, 
                                        non_deg_genes = non_deg_genes,
                                        genes = genes,
                                        clique_significance)
    }
    
    significant_vector <- apply(X = size_and_deg, MARGIN = 2, FUN = check_significance, cutoff_per_size = cutoff_per_size)
    
    significant_cliques <- unique(c(significant_cliques, unlist(clique_genes[significant_vector])))
  }
  
  
  significant_cliques <- unname(unique(unlist(sapply(X = significant_cliques, 
                                                     FUN = index_to_entrez, 
                                                     name_table = name_table))))
  
  close(con)
  
  file.remove(clique_file)
  
  return (unique(significant_cliques))
  
}

entrez_to_index <- function(entrez_id, name_table){
  name_table[which(name_table[,1] == entrez_id), 2 ]
}
index_to_entrez <- function(gene_index, name_table){
  name_table[which(name_table[,2] == gene_index), 1 ]
}

get_clique_cutoffs <- function(clique_size, deg_names, non_deg_genes, genes, clique_significance){
  
  create_fisher_tables(clique_size = clique_size, 
                       total_n_genes = length(genes), 
                       n_sig = length(deg_names), 
                       unsig = length(non_deg_genes)) %>% 
    evaluate_table(fisher_table = ., 
                   clique_significance = clique_significance, 
                   min_deg_in_clique = 1)
}

check_significance <- function(clique, cutoff_per_size){
  if (clique[2] >= cutoff_per_size[clique[1]]){
    return (TRUE)
  }else{
    return(FALSE)
  }
}

scale_observations <- function(sf, scores){
  scores * (length(scores) * sf) / sum(scores)
}
