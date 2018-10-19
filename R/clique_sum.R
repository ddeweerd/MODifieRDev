#' Clique Sum
#'
#' An implementation of the clique-based disease module inference method proposed by Barrenäs et al. 
#'
#' @param MODifieR_input A MODifieR input object produced by \code{\link{create_input}} function
#' @param ppi_network A PPi network as a 2 column dataframe
#' @param n_iterations Number of iterations to be performed for the permutation based p-value
#' @param clique_significance p-value for cliques to be considered significant
#' @param deg_cutoff  p-value cutoff for differentialy expressed genes
#' @param min_clique_size Minimal size for cliques 
#' @param min_deg_in_clique Minimum number of DEGs to be present in the clique
#' @return clique_sum returns an object of class "MODifieR_module" with subclass "Clique_Sum". 
#' This object is a named list containing the following components:
#' \item{module_genes}{A character vector containing the genes in the final module}
#' \item{settings}{A named list containing the parameters used in generating the object}
#' @details 
#' The clique sum algorithm obtains \emph{maximal cliques}, cliques that are not contained within another clique,
#'  from the \code{ppi_network}. The minimal size of the obtained cliques is equal to \code{min_clique_size}. 
#'  The union of maximal cliques with a permutation based p-value below \code{clique_significance} and at 
#'  least \code{min_deg_in_clique} is the disease module
#' @references 
#' \cite{Barrenäs F, Chavali S, Alves AC, et al. Highly interconnected genes in disease-specific networks are enriched for disease-associated polymorphisms. Genome Biology. 2012;13(6):R46. doi:10.1186/gb-2012-13-6-r46.}
#' @export
clique_sum <- function(MODifieR_input, ppi_network, simplify_graph = T, n_iterations = 10000, clique_significance = 0.05,
                       deg_cutoff = 0.05, min_clique_size = 5, min_deg_in_clique = 3, dataset_name = NULL){
  # Retrieve settings
  default_args <- formals()
  user_args <- as.list(match.call(expand.dots = T)[-1])
  settings <- c(user_args, default_args[!names(default_args) %in% names(user_args)])
  
  if (!is.null(dataset_name)){
    settings$MODifieR_input <- dataset_name
  }
  
  deg_genes <- MODifieR_input$diff_genes
  
  #Convert 2 column dataframe to named vector, names will be the gene names, values will be the p value
  deg_genes <- dataframe_to_vector(as.data.frame(deg_genes))
  
  deg_genes <- deg_genes[deg_genes < cutoff]
  #Convert PPi to igraph object
  graphed_frame <- igraph::graph.data.frame(unique(ppi_network) , directed = FALSE)
  #Only retain DEGs that are also present in the PPi network
  deg_genes <- deg_genes[intersect(names(deg_genes), igraph::V(graphed_frame)$name)]

  if (simplify_graph == T){
    graphed_frame <- igraph::simplify(graphed_frame)
  }
  #Get maximal cliques from the PPi network
  cv <- igraph::max_cliques(graph = graphed_frame, min = min_clique_size)
  #Get number of DEGs per clique
  n_deg_clique <- sapply(X = cv, FUN = function(x){length(which(names(x)%in% names(deg_genes)))})
  #Retrieve the indici of cliques which have at least min_deg_in_clique number of DEGs in them
  min_deg <- n_deg_clique >= min_deg_in_clique
  #Subset the clique list and the vector with number of DEGs per clique
  cv <- cv[min_deg]
  n_deg_clique <- n_deg_clique[min_deg]
  #Sum the p values for each gene for each clique. Will return a vector with summed p-values for each clique
  unpermuted_scores <- sapply(X=cv ,FUN = function(x){sum(deg_genes[names(x)], na.rm = T)} ,simplify = TRUE)
  #Sort the number of DEGs in each clique and only take the unique values 
  deg_sizes <- sort(unique(n_deg_clique))
  #Get p values for each clique
  p_values <- unlist(sapply(X = deg_sizes, FUN = permute_scores,
                            n_iterations = n_iterations, deg_genes = deg_genes, unpermuted_scores,
                            n_deg_clique = n_deg_clique))
  #Retain cliques that have a significance < than the cutoff
  significant_cliques <- p_values < cutoff

  cv <- cv[significant_cliques]
  #Unlist all cliques and remove duplicate genes. This is the final disease module
  final_genenames <-unique(names(unlist(cv)))
  #Construct list and give it the correct class
  new_clique_sum_module <- list("module_genes" =  final_genenames,
                                "settings" = settings)

  class(new_clique_sum_module) <- c("MODifieR_module", "Clique_Sum")

  return(new_clique_sum_module)
}
#Take n_iterations times randomly a number DEGs equal to the size of deg_size, and sum these p-values
#compare the random summed p values to the true summed p value for that clique and divide by n_iterations to
#get the p-value
permute_scores <- function(deg_size, n_iterations, deg_genes, unpermuted_scores, n_deg_clique){
  permutations <- replicate(n = n_iterations, expr = sum(deg_genes[sample(length(deg_genes), deg_size, replace = F)]))
  p_values <- sapply(unpermuted_scores[n_deg_clique == deg_size], FUN = score_p_value,
                     permuted_scores = permutations, n_iterations = n_iterations)
  return(p_values)
}
#
score_p_value <- function(x, permuted_scores, n_iterations){
  (sum(permuted_scores > x)) / n_iterations
}
