#' DiffCoEx
#' 
#' An implementation of the DiffCoEx co-expression based algorithm
#' @inheritParams clique_sum
#' @inheritParams flashClust::flashClust
#' @inheritParams dynamicTreeCut::cutreeDynamic
#' @inheritParams WGCNA::mergeCloseModules
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

diffcoex <- function(MODifieR_input, beta = 6, cor_method = "spearman",
                     cluster_method = "average", cuttree_method = "hybrid",
                     cut_height = 0.996, deepSplit = 0, pamRespectsDendro = F,
                     minClusterSize = 20, cutHeight = 0.2, dataset_name = NULL){
  
  # Retrieve settings
  default_args <- formals()
  user_args <- as.list(match.call(expand.dots = T)[-1])
  settings <- c(user_args, default_args[!names(default_args) %in% names(user_args)])
  
  if (!is.null(dataset_name)){
    settings$MODifieR_input <- dataset_name
  }
  #Get relevant input data from input object 
  dataset1 <- t(MODifieR_input$annotated_exprs_matrix[,MODifieR_input$group_indici[[1]]])
  dataset2 <- t(MODifieR_input$annotated_exprs_matrix[,MODifieR_input$group_indici[[2]]])

  AdjMatC1<-sign(stats::cor(dataset1, method = cor_method))*(stats::cor(dataset1, method = cor_method))^2
  AdjMatC2<-sign(stats::cor(dataset2, method = cor_method))*(stats::cor(dataset2, method = cor_method))^2
  diag(AdjMatC1)<-0
  diag(AdjMatC2)<-0
  
  WGCNA::collectGarbage()
  
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
  #Build new diffxoex object
  new_diffcoex_module <- construct_diffcoex_module(module_list = module_genes_list,
                                                   annotation_table = MODifieR_input$annotation_table,
                                                   color_vector = color_vector,
                                                   settings = settings)
  
  return (new_diffcoex_module)
  
}
#Function to aggregate colors in a list
aggregate_colors <- function(color, genes, color_vector){
  genes[which(names(color_vector) == color)]
}
#Module constructor function
construct_diffcoex_module <- function(module_list, annotation_table, color_vector, settings){
  
  module_genes <- as.vector(unlist(unname(module_list)))
  module_colors <- names(module_list)
  
  new_diffcoex_module <- list("module_genes" =  module_genes,
                              "module_colors"= module_colors,
                              "color_vector" =  color_vector,
                              "settings" = settings)
  
  class(new_diffcoex_module) <- c("MODifieR_module", "DiffCoEx")
  
  
  return (new_diffcoex_module)
}
