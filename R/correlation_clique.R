#' Clique_correlation
#'
#' A clique based method to find a disease module from correlated gene expression
#'
#' @inheritParams clique_sum
#' @param iteration Number of iterations to be performed
#' @param probabilityScaleFactor Scale for enriched cliques 
#' @param frequency_cutoff significance out of number of iterations performed (default = 0.05) 
#' @param signif_cutoff Cutoff for Fisher exact test for cliques
#' @param deg_cutoff Cutoff for significant genes
#' @details 
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
#' @return correlation_clique returns an object of class "MODifieR_module" with subclass "Correlation_clique". 
#' This object is a named list containing the following components:
#' \item{module_genes}{A character vector containing the genes in the final module}
#' \item{frequency_table}{A table containing the fraction of times the genes were present in an iteration module}
#' \item{settings}{A named list containing the parameters used in generating the object}
#' @export
correlation_clique <- function(MODifieR_input, ppi_network,
                                frequency_cutoff = .5, cutOffRank = 700,
                                probabilityScaleFactor = 0.6,
                                iteration = 50, signif_cutoff = 0.01, 
                                deg_cutoff = 0.05,
                                dataset_name = NULL){
  
  # Retrieve settings
  default_args <- formals()
  user_args <- as.list(match.call(expand.dots = T)[-1])
  settings <- c(user_args, default_args[!names(default_args) %in% names(user_args)])
  
  if (!is.null(dataset_name)){
    settings$MODifieR_input <- dataset_name
  }
  
  springConnection <- ppi_network[,1:3]
  
  pValueMatrix <- MODifieR_input$diff_genes
  
  pValueMatrix <- na.omit(pValueMatrix)

  overlap <- springConnection[,1] %in% pValueMatrix[,1]
  
  springConnection <- springConnection[overlap,]
  
  overlap <- springConnection[,2] %in% pValueMatrix[,1]
  
  springConnection <- springConnection[overlap,]
  
  overlap <- springConnection[,1] %in% rownames(MODifieR_input$annotated_exprs_matrix)
  
  springConnection <- springConnection[overlap,]
  
  overlap <- springConnection[,2] %in% rownames(MODifieR_input$annotated_exprs_matrix)
  
  springConnection <- springConnection[overlap,]
  
  springConnection <- springConnection[springConnection[,3] > cutOffRank,]
  
  springConnection[,3] <- springConnection[,3] / 1000
  #Calculates (1- p value) for all relevant connections in the PPi network, and is stored in column 1.
  #Column 2 will contain the PPi score from the PPi network
  corrPvalues <- cbind(apply(X = springConnection, MARGIN = 1, FUN = calculate_correlation, 
                             expression_matrix = MODifieR_input$annotated_exprs_matrix),
                       springConnection[,3])
  #Calculates square root of p-value * PPi score. This is then scaled by probabilityScaleFactor. 
  pval_score <- apply(X = corrPvalues, MARGIN = 1, FUN = function(x){sqrt(x[1] * x[2]) * probabilityScaleFactor})
  #Get the significant genes
  signifgenes <- pValueMatrix$gene[pValueMatrix$pvalue < deg_cutoff]
  n_unsignif <- nrow(pValueMatrix)
  

  module_list <- list() 
  #For every iteration:
  for (i in 1:iteration){
    #Compare computed scores to random scores 
    to_pass <- pval_score > runif(n = length(pval_score)) * 1
    #Get a graphed object for the edges that have a higher score than random (to_pass vector)
    graphed_adjeceny <- graph_score(adjecency_list = cbind(springConnection[,1:2], to_pass))
    #Get maximal cliques from the graphed object
    cliques <- igraph::max_cliques(graph = graphed_adjeceny, min = 2)
    #Infer modules using the maximal cliques
    module_list[[i]] <- infer_module(cliques = cliques, signifgenes = signifgenes, n_unsignif = n_unsignif, signif_cutoff = signif_cutoff )
  }
  #Get the frequency of each gene that has been present in a module, expressed in fraction
  tabled_frequencies <- table(unlist(module_list)) / iteration
  #The resulting module is going to consist of all genes that have are seen more frequently in the iterations
  #than the frequency cutoff (default 0.5, so present in 50% of all iterations) 
  module_genes <- names(tabled_frequencies[tabled_frequencies >= frequency_cutoff])
  
  new_correlation_clique_module <- list("module_genes" = module_genes,
                                        "frequency_table" = tabled_frequencies,
                                        "settings" = settings)
  
  class( new_correlation_clique_module) <- c("MODifieR_module", "Correlation_clique")
  
  return( new_correlation_clique_module)
}

#Calculate fisher exact test for each clique. Returns a vector with all 
#genes in significant clique
infer_module <- function(cliques, signifgenes, n_unsignif, signif_cutoff){
  #Total number of significant genes
  n_signif <- length(signifgenes)
  #Total number of not significant genes
  n_unsignif <- n_unsignif - n_signif
  #Clique table, matrix with dimensions: 4 columns, number of cliques = number of rows:
  #Column 1 Number of significant genes in clique
  #Column 2 Number of unsignificant genes in clique
  #Column 3 Number of significant genes NOT in clique
  #Column 4 Number of unsignficant genes NOT in clique
  
  #Columns 1+2 
  clique_table <- cbind(sapply(X = cliques, FUN = function(x) sum(as.vector(x) %in% signifgenes)), 
                        sapply(X = cliques, FUN = function(x) sum(!as.vector(x) %in% signifgenes)))
  #cbind column 3
  clique_table <- cbind(clique_table, n_signif - clique_table[,1])
  #Get indici of cliques with at least 1 significant gene
  #Avoids calculating Fisher exact for cliques with no significant gene
  cliques_with_signif <- clique_table[ ,1] != 0

  #cbind column 4
  clique_table <- cbind(clique_table, n_unsignif - clique_table[,2])
  
  #Subset clique table w
  clique_table <- clique_table[cliques_with_signif,]
  #Subset cliques
  cliques <- cliques[cliques_with_signif]
  
  #Get fisher exact test p values for each cliqiue. A, B, C and D are in column 1:4
  clique_pvals <- apply(X = clique_table, MARGIN = 1, FUN = fisher_pval)
  
  #Logical vector, TRUE is significant
  indici_signif_cliques <- clique_pvals < signif_cutoff
  #Unlist, unique and vectorize all significant cliques
  module <- unique(names(unlist(cliques[indici_signif_cliques])))
  
  return (module)
  
}
#Pval Helper function
fisher_pval <- function(clique_row){
  fisher.test(x = matrix(data = clique_row, nrow = 2, byrow = T), alternative = "g")$p.value
}
#This function takes as the first argument a row from a PPi network, retrieves the corresponding row and column
#by row and column name and computes (1- correlation P value). 
calculate_correlation <- function(row, expression_matrix){
  row <- as.character(row)
  
  x_y <- c(which(rownames(expression_matrix) == row[1]), which(rownames(expression_matrix) == row[2]))
  
  1 - round(cor.test(x = expression_matrix[x_y[1],], y = expression_matrix[x_y[2],])$p.value, 3)
}
#Convert edge list to igraph object, and convert to adjacency matrix. This matrix will be converted to an adjacency graph.
#As this adjacency matrix will be used to search for maximal cliques, the mode is undirected
graph_score <- function(adjecency_list){
  g <- igraph::graph.data.frame(adjecency_list)
  score_matrix <- as.matrix(igraph::get.adjacency(g , attr = colnames(adjecency_list)[3]))
  
  #Remove rows that do not have a score of at least 1 in scoreMatrix
  rows_to_keep <- rowSums(x = score_matrix) != 0
  
  score_matrix <- score_matrix[rows_to_keep, rows_to_keep]
  #Set lower triangle to 0. 
  score_matrix[lower.tri(score_matrix)] <- 0
  
  z <- igraph::graph.adjacency(adjmatrix = score_matrix, mode = "undirected")
  z1 <- igraph::simplify(z)
  
  return(z1)
}

