#' Clique_correlation
#'
#' A clique based method to find a disease module from correlated gene expression
#'
#'@param pval_matrix A nx2 dataframe of the entrez ids and their respective p-values
#'@param correlation_matrix A nxn matrix of Correlated p-values of n genes 
#'@param ppi_network A dataframe of PPI network of your choice
#'@param iteration Number of iterations to be performed
#'@param cutoffRank Score cutoff of the string PPI network
#'@param cliquetestCriteria type of test to be perfomed for enrichment of cliques (default = "FisherTest")
#'@param probabailitySacleFactor Scale for enriched cliques 
#'@param frequency_cutoff Numeric significance out of number of iterations performed (default = 0.05) 
#'@param signif_cutoff Cutoff for Fisher exact test for cliques
#'@param gene_significance Cutoff for significant genes
#'@return A MODifieR class object with disease module and settings
#'@export
correlation_clique <- function(MODifieR_input, ppi_network,
                                frequency_cutoff = .5, cutOffRank = 700,
                                probabilityScaleFactor = 0.6,
                                iteration = 50, signif_cutoff = 0.01, 
                                gene_significance = 0.05,
                                dataset_name = NULL){
  
  # Retrieve settings
  default_args <- formals()
  user_args <- as.list(match.call(expand.dots = T)[-1])
  settings <- c(user_args, default_args[!names(default_args) %in% names(user_args)])
  
  if (!is.null(dataset_name)){
    settings$MODifieR_input <- dataset_name
  }
  
  springConnection <- ppi_network[,1:3]
  
  pValueMatrix <- plyr::ddply(.data = MODifieR_input$diff_genes, 
                           .variables = "ENTREZID", .fun = plyr::summarise, pvalue = min(P.Value))
  
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
  #Calculates (1- p value) for all relevent connections in the PPi network, and is stored in column 1.
  #Column 2 will contain the PPi score from the PPi network
  corrPvalues <- cbind(apply(X = springConnection, MARGIN = 1, FUN = calculate_correlation, 
                             expression_matrix = MODifieR_input$annotated_exprs_matrix),
                       springConnection[,3])
  #Calculates square root of P value * PPi score. This is then scaled by probabilityScaleFactor. 
  pval_score <- apply(X = corrPvalues, MARGIN = 1, FUN = function(x){sqrt(x[1] * x[2]) * probabilityScaleFactor})
  #Get the significant genes
  signifgenes <- pValueMatrix$gene[pValueMatrix$pvalue < gene_significance]
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
                                        "settings" = settings,
                                        "frequency_table" = tabled_frequencies)
  
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
#As this adjecency matrix will be used to search for maximal cliques, the mode is undirected
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

