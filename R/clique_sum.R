#' Clique Sum
#'
#' A Clique based Algorithm proposed by Fredrik Barranais to produce a disease module from Differentially expressed genes using a core sum approach
#'
#' @param deg_genes A nx2 dataframe consisting of differentially expressed genes and their respective -log10p values
#' @param igraph_network A path to dataframe csv PPi network of your choice
#' @param n_iter number of iterations to be performed
#' @param cutoff  Desired P-value cutoff
#' @param min_clique_size Numeric shortest clique size (default = 5)
#' @param min_deg_clique Numeric Minimum degree of the clique
#' @return A MODifieR class object with disease module and settings
#' @export
clique_sum <- function(MODifieR_input, ppi_network, simplify_graph = T, n_iter = 10000,
                       deg_cutoff = 0.05, min_clique_size = 5, min_deg_clique = 3, dataset_name = NULL){
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
  
  graphed_frame <- igraph::graph.data.frame(unique(ppi_network) , directed = FALSE)

  deg_genes <- deg_genes[intersect(names(deg_genes), V(graphed_frame)$name)]

  if (simplify_graph == T){
    graphed_frame <- igraph::simplify(graphed_frame)
  }

  cv <- igraph::max_cliques(graph = graphed_frame, min = min_clique_size, subset = names(deg_genes))

  n_deg_clique <- sapply(X = cv, FUN = function(x){length(which(names(x)%in% names(deg_genes)))})
  min_deg <- n_deg_clique >= min_deg_clique

  cv <- cv[min_deg]
  n_deg_clique <-n_deg_clique[min_deg]

  unpermuted_scores <- sapply(X=cv ,FUN = function(x){sum(deg_genes[names(x)], na.rm = T)} ,simplify = TRUE)
  deg_sizes <- sort(unique(n_deg_clique))

  p_values <- unlist(sapply(X = deg_sizes, FUN = permute_scores,
                            n_iter = n_iter, deg_genes = deg_genes, unpermuted_scores,
                            n_deg_clique = n_deg_clique))

  significant_cliques <- p_values < cutoff

  cv <- cv[significant_cliques]

  final_genenames <-unique(names(unlist(cv)))

  new_clique_sum_module <- list("module_genes" =  final_genenames,
                                "settings" = settings)

  class(new_clique_sum_module) <- c("MODifieR_module", "Clique_Sum")

  #new_clique_sum_module <- symbol_to_entrez(new_clique_sum_module)

  return(new_clique_sum_module)
}

permute_scores <- function(deg_size, n_iter, deg_genes, unpermuted_scores, n_deg_clique){
  permutations <- replicate(n = n_iter, expr = sum(deg_genes[sample(length(deg_genes), deg_size, replace = F)]))
  p_values <- sapply(unpermuted_scores[n_deg_clique == deg_size], FUN = score_p_value,
                     permuted_scores = permutations, n_iter = n_iter)
  return(p_values)
}

score_p_value <- function(x, permuted_scores, n_iter){
  (sum(permuted_scores > x)) / n_iter
}
