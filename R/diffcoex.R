#' DiffCoEx
#' 
#' An implementation of the DiffCoEx co-expression based algorithm
#' @inheritParams clique_sum_exact
#' @inheritParams flashClust::flashClust
#' @inheritParams dynamicTreeCut::cutreeDynamic
#' @inheritParams WGCNA::mergeCloseModules
#' @inheritParams wgcna
#' @param cut_height Maximum joining heights that will be considered. 
#' @param cluster_method the agglomeration method to be used. 
#' This should be (an unambiguous abbreviation of) one of "ward", 
#' "single", "complete", "average", "mcquitty", "median" or "centroid".
#' This applies to hierachical clustering.
#' @param beta User-defined soft thresholding power
#' For method=="tree" it defaults to 0.99. For method=="hybrid" it 
#' defaults to 99% of the range between the 5th percentile and the
#'  maximum of the joining heights on the dendrogram.
#' @param cor_method a character string indicating which correlation coefficient (or covariance)
#' is to be computed. One of "pearson" (default), "kendall", or "spearman": can be abbreviated.
#' @param cuttree_method Chooses the method to use. Recognized values are "hybrid" and "tree". 
#' 
#' @details
#'  DiffCoEx is a method for identifying correlation pattern changes, which builds on the commonly 
#'  used Weighted Gene Coexpression Network Analysis (WGCNA) framework for coexpression analysis.
#' 
#' @return diffcoex returns an object of class "MODifieR_module" with subclass "DiffCoEx". 
#' This object is a named list containing the following components:
#' \item{module_genes}{A character vector containing the genes in the final module}
#' \item{module_colors}{A character vector containing the colors that make up the final disease module}
#' \item{color_vector}{A named character vector containing the genes as values and the color as name}
#' \item{settings}{A named list containing the parameters used in generating the object}
#' @references 
#' \cite{Tesson, B. M., Breitling, R., & Jansen, R. C. (2010). DiffCoEx: a simple and sensitive
#'  method to find differentially coexpressed gene modules. BMC Bioinformatics, 11, 497.
#'   \url{https://doi.org/10.1186/1471-2105-11-497}}
#' 
#' @export

diffcoex <- function(MODifieR_input, beta = NULL, cor_method = "spearman",
                     cluster_method = "average", cuttree_method = "hybrid",
                     cut_height = 0.996, deepSplit = 0, pamRespectsDendro = F,
                     minClusterSize = 20, cutHeight = 0.2, 
                     pval_cutoff = 0.05, dataset_name = NULL, deg_cutoff = 0.05){
  
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
  if (is.null(beta)){
    beta <- get_softthreshold_value(datExpr = t(MODifieR_input$annotated_exprs_matrix), 
                                    nSamples = ncol(MODifieR_input$annotated_exprs_matrix))
    settings$beta <- beta
  }
  #Get relevant input data from input object 
  dataset1 <- t(MODifieR_input$annotated_exprs_matrix[,MODifieR_input$group_indici[[1]]])
  dataset2 <- t(MODifieR_input$annotated_exprs_matrix[,MODifieR_input$group_indici[[2]]])
  
  AdjMatC1<-sign(stats::cor(dataset1, method = cor_method))*(stats::cor(dataset1, method = cor_method))^2
  AdjMatC2<-sign(stats::cor(dataset2, method = cor_method))*(stats::cor(dataset2, method = cor_method))^2
  diag(AdjMatC1)<-0
  diag(AdjMatC2)<-0
  
  dissTOMC1C2 <- WGCNA::TOMdist((abs(AdjMatC1-AdjMatC2)/2)^(beta/2))
  WGCNA::collectGarbage()
  
  #Hierarchical clustering is performed using the Topological Overlap of the adjacency difference as input distance matrix
  geneTreeC1C2 <- flashClust::flashClust(stats::as.dist(dissTOMC1C2), method = cluster_method);
  
  #We now extract modules from the hierarchical tree. This is done using cutreeDynamic. Please refer to WGCNA package documentation for details
  dynamicModsHybridC1C2 <- dynamicTreeCut::cutreeDynamic(dendro = geneTreeC1C2,
                                                         distM = dissTOMC1C2,
                                                         method = cuttree_method,
                                                         cutHeight = cut_height,
                                                         deepSplit = deepSplit,
                                                         pamRespectsDendro = pamRespectsDendro,
                                                         minClusterSize = minClusterSize);
  
  #Every module is assigned a color. Note that GREY is reserved for genes which do not belong to any differentially coexpressed module
  dynamicColorsHybridC1C2 <- WGCNA::labels2colors(dynamicModsHybridC1C2)
  
  #the next step merges clusters which are close (see WGCNA package documentation)
  colorh1C1C2 <- WGCNA::mergeCloseModules(rbind(dataset1,dataset2),
                                          dynamicColorsHybridC1C2,
                                          cutHeight=cutHeight)$color
  
  color_vector <- rownames(MODifieR_input$annotated_exprs_matrix)
  
  names(color_vector) <- colorh1C1C2
  #Retrieve colors
  colors <- unique(colorh1C1C2)
  #But exclude grey
  colors <- colors[colors!="grey"]
  #Agrregate all module color genes in a list
  module_genes_list <- lapply(colors, FUN = aggregate_colors, 
                              genes = rownames(MODifieR_input$annotated_exprs_matrix), 
                              color_vector = color_vector)
  #Name this list according to color
  names(module_genes_list) <- colors
  
  pval_null_distr <- lapply(X = colors, FUN = diffcoex_get_p_values, dataset1 = dataset1, 
                            dataset2 = dataset2, 
                            group1_indici = MODifieR_input$group_indici[[1]], 
                            color_vector = color_vector, cor_method = cor_method)
  
  module_p_values <- cbind(sapply(X = pval_null_distr, function(x){x$emp_pvalue}),
                           sapply(X = pval_null_distr, function(x){x$z_score}))
  
  null_distributions <- t(sapply(X = pval_null_distr, function(x)x$null_distribution))
  
  colnames(module_p_values) <- c("p_values", "z_scores")
  
  rownames(null_distributions) <- colors
  
  fisher_exact_p_values <- lapply(X = module_genes_list, FUN = color_module_exact_test, 
                                  MODifieR_input = MODifieR_input, deg_cutoff = deg_cutoff)
  
  fisher_exact_p_values <- do.call(rbind, fisher_exact_p_values)
  
  module_p_values <- cbind(module_p_values, fisher_exact_p_values)
  module_genes <- as.vector(unlist(unname(module_genes_list[which(module_p_values < pval_cutoff)])))
  module_colors <- names(module_genes_list[which(module_p_values < pval_cutoff)])
  
  new_diffcoex_module <- construct_diffcoex_module(module_genes = module_genes,
                                                   module_colors = module_colors,
                                                   module_p_values = module_p_values,
                                                   color_vector = color_vector,
                                                   null_distributions = null_distributions,
                                                   settings = settings)
  
  return (new_diffcoex_module)
  
}
#Function to aggregate colors in a list
aggregate_colors <- function(color, genes, color_vector){
  genes[which(names(color_vector) == color)]
}
#Module constructor function
construct_diffcoex_module <- function(module_genes, module_colors, module_p_values,
                                      color_vector, null_distributions, settings){
  
  
  new_diffcoex_module <- list("module_genes" =  module_genes,
                              "module_colors" = module_colors,
                              "module_p_values" = module_p_values,
                              "color_vector" =  color_vector,
                              "null_distributions" = null_distributions,
                              "settings" = settings)
  
  class(new_diffcoex_module) <- c("MODifieR_module", "DiffCoEx")
  
  
  return (new_diffcoex_module)
}
#' Returns new DiffCoEx module objects by color
#' @param diffcoex_module Module object that has been produced by \code{diffcoex} function
#' @details 
#' The  \code{DiffCoEx} module object is split into a series of \code{DiffCoEx} objects by color.
#' Eevery significant color in the module will be its own \code{DiffCoEx} module object
#' 
#' @return 
#' 
#' A list of \code{DiffCoEx} module objects
#' 
#' @seealso
#' 
#' \code{\link{diffcoex}}
#' 
#' @export
diffcoex_split_module_by_color <- function(diffcoex_module){
  module_list <- list()
  color_vector <- diffcoex_module$color_vector
  for (color in diffcoex_module$module_colors){
    module_genes <- unname(color_vector[which(names(color_vector) == color)])
    module_list[[color]] <- construct_diffcoex_module(module_genes = module_genes,
                                                      module_colors = color,
                                                      color_vector = color_vector,
                                                      module_p_values = diffcoex_module$module_p_values,
                                                      settings = diffcoex_module$settings)
  }
  return(module_list)
}

diffcoex_get_p_values <- function(color, dataset1, dataset2, group1_indici, color_vector, cor_method){
  
  scaled_combined_dataset <- rbind(scale(dataset1),scale(dataset2))
  
  pval_null_distr <- permute_modules(color1 = color, color2 = color, 
                                     dataset = scaled_combined_dataset, n_samples = length(group1_indici), 
                                     color_vector = color_vector, sample_labels = group1_indici, cor_method = cor_method)
  
  return(pval_null_distr)
}

permute_modules <- function(color1, color2, dataset, n_samples, color_vector, sample_labels, cor_method){
  
  module_score <- module_dispersion(color1 = color1, color2 = color2, dataset = dataset, 
                                    n_samples = n_samples, color_vector = color_vector, 
                                    sample_labels = sample_labels, cor_method = cor_method)
  
  null_distribution <- replicate(n = 1000, expr = module_dispersion(color1 = color1, 
                                                                    color2 = color2, 
                                                                    dataset = dataset, 
                                                                    n_samples = n_samples, 
                                                                    color_vector = color_vector,
                                                                    cor_method = cor_method))
  
  emp_pvalue <- sum(null_distribution >= module_score) / 1000
  
  z_score <- (module_score - mean(null_distribution)) / sd(null_distribution)
  
  return(list("emp_pvalue" = emp_pvalue, "z_score" = z_score, "null_distribution" =  null_distribution))
}

module_dispersion <- function(color1, color2, dataset, n_samples, color_vector, sample_labels = NULL, cor_method){
  if (is.null(sample_labels)){
    sample_labels <- sample(size = n_samples, x = nrow(dataset))
  }
  
  if (color1 == color2){
    
    module1 <- dataset[sample_labels, unname(color_vector[names(color_vector) == color1])]
    module2 <- dataset[-sample_labels, unname(color_vector[names(color_vector) == color2])]
    
    difCor <- (cor(module1, method = cor_method) - cor(module2, method = cor_method)) ^2
    n <- ncol(module1)
    return  ((1/((n^2 -n)/2) * (sum(difCor)/2))^(.5))
  }else{
    module1 <- dataset[sample_labels, unname(color_vector[names(color_vector) == color1])]
    module2 <- dataset[-sample_labels, unname(color_vector[names(color_vector) == color1])]
    module3 <- dataset[sample_labels, unname(color_vector[names(color_vector) == color2])]
    module4 <- dataset[-sample_labels, unname(color_vector[names(color_vector) == color2])]
    
    difCor <-(cor(module1,module3, method = cor_method) - cor(module2, module4, method = cor_method))^2
    n1 <- length(unname(color_vector[names(color_vector) == color1]))
    n2 <- length(unname(color_vector[names(color_vector) == color2]))
    
    return((1/((n1*n2)) * (sum(difCor)))^(.5))
    
  }
  
}

color_module_exact_test <- function(color_module, MODifieR_input, deg_cutoff){
  all_genes <- intersect(MODifieR_input$diff_genes$gene, 
                         rownames(MODifieR_input$annotated_exprs_matrix))
  diff_genes <- MODifieR_input$diff_genes[MODifieR_input$diff_genes$gene %in% all_genes, ]
  
  deg_genes <- diff_genes$gene[diff_genes$pvalue < deg_cutoff]
  non_deg_genes <- diff_genes$gene[diff_genes$pvalue >= deg_cutoff]
  
  color_module <- intersect(all_genes, color_module)
  co_expression_enrichment <- c(
    sum(color_module %in% deg_genes),
    sum(color_module %in% non_deg_genes),
    length(deg_genes) - sum(color_module %in% deg_genes),
    length(non_deg_genes) - sum(color_module %in% non_deg_genes))
  
  expected_degs <- (length(deg_genes) / length(all_genes)) * length(color_module)
  observed_degs <- sum(color_module %in% deg_genes)
  
  return(c("Fisher_pval" = module_exact_test <- MODifieRDev:::fisher_pval(co_expression_enrichment), 
           "expected_DEGs" = expected_degs, "Observed_DEGs" = observed_degs))
}
