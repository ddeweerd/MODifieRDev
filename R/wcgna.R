#Module object constructor
wgcna_module_constructor <- function(module_genes, info_table, 
                                     module_cor_and_p_value, powerEstimate,
                                     module_colors, input_class, settings){
  
  new_wgcna_module <- list("module_genes" =  module_genes,
                           "info_table" = info_table,
                           "correlation_to_trait_table" = module_cor_and_p_value,
                           "softthreshold_value" = powerEstimate,
                           "module_colors" = module_colors,
                           "settings" = settings)
  
  class(new_wgcna_module) <- c("MODifieR_module", "WGCNA", input_class)
  
  return(new_wgcna_module)
}

#' An implementation of WGCNA to correlate coexpression modules to disease
#'
#' @inheritParams clique_sum_exact
#' @inheritParams WGCNA::blockwiseModules
#' @inheritParams moda
#' @param pval_cutoff The p-value cutoff to be used for significant co-expression modules (colors)
#' 
#' @details
#' wgcna is an implementation of WGCNA that
#' associates co-expression modules (denoted by color) to a trait. Co-expression modules with an
#' adjusted p-value < \code{pval_cutoff} will make up the final disease module.
#' 
#' The algorithm infers co-expression modules from combined expression dataset from both \code{group1}
#' and \code{group2}.
#' Co-expression modules are then correlated to trait (group 1 ~ group 2). 
#' 
#' After analysis there are some post-processing functions available:
#' 
#' \itemize{
#' \item{\code{\link{wgcna_get_all_module_genes}}} Get a list with all genes sorted by module color
#' \item{\code{\link{wgcna_get_module_genes_by_sign}}} Get a module with either only postively correlated
#' genes or negatively correlated genes
#' \item{\code{\link{wgcna_adjust_significance}}} Adjust p-value cutoff
#' \item{\code{\link{wgcna_split_module_by_color}}} Get a list where each color is a separate module
#' \item{\code{\link{wgcna_set_module_size}}} Get a module close to a specific size
#' }
#' 
#' @return wgcna returns an object of class "MODifieR_module" with subclass "WGCNA". 
#' This object is a named list containing the following components:
#' \item{module_genes}{A character vector containing the genes in the final module}
#' \item{info_table}{A data.frame containing all genes and their assigned colors}
#' \item{correlation_to_trait_table}{A data.frame containing all module colors and their p- and adjusted p-value}
#' \item{softthreshold_value}{A numeric, the soft threshold power that is used. See: \code{\link[WGCNA]{pickSoftThreshold}}}
#' \item{module_colors}{A character vector containing the colors that make up the final disease module}
#' \item{settings}{A named list containing the parameters used in generating the object}
#' 
#' @references 
#' \cite{Langfelder P and Horvath S, WGCNA: an R package for weighted correlation network analysis. 
#' BMC Bioinformatics 2008, 9:559 \url{doi:10.1186/1471-2105-9-559} }
#' 
#' 
#' \cite{Peter Langfelder, Steve Horvath (2012). Fast R Functions for Robust Correlations and Hierarchical Clustering. J
#' ournal of Statistical Software, 46(11), 1-17. URL \url{http://www.jstatsoft.org/v46/i11/}}
#' 
#' @author Dirk de Weerd
#' 
#' @export
wgcna <- function(MODifieR_input, group_of_interest, minModuleSize = 30, deepSplit = 2, pamRespectsDendro = F,
                  mergeCutHeight = 0.1, numericLabels = T,  pval_cutoff = 0.05, corType = "bicor",
                  maxBlockSize = 5000, TOMType = "signed", saveTOMs = T, maxPOutliers = 0.1, 
                  dataset_name = NULL){
  
  # Retrieve settings
  evaluated_args <- c(as.list(environment()))
  settings <- as.list(stackoverflow::match.call.defaults()[-1])
  replace_args <- names(settings)[!names(settings) %in% unevaluated_args]
  for (argument in replace_args) {
    settings[[which(names(settings) == argument)]] <- evaluated_args[[which(names(evaluated_args) == 
                                                                              argument)]]
  }
  #Validate the input parameters
  check_expression_matrix(MODifieR_input)
  validate_inputs(settings)
  if (!(group_of_interest == 1 | group_of_interest == 2)){
    stop("Group of interest is not 1 or 2", call. = F)
  }
  
  if (!is.null(dataset_name)){
    settings$MODifieR_input <- dataset_name
  }
  
  allowWGCNAThreads()
  
  traits <- wgcna_get_trait_data(MODifieR_input = MODifieR_input, group_of_interest = group_of_interest)
  
  datExpr <- t(MODifieR_input$annotated_exprs_matrix)
  
  nGenes <- ncol(datExpr)
  nSamples <- nrow(datExpr)
  
  powerEstimate <- get_softthreshold_value(datExpr = datExpr, nSamples = nSamples)
  
  net <- blockwiseModules(datExpr, power = powerEstimate,
                          TOMType = TOMType, minModuleSize = minModuleSize,
                          mergeCutHeight = mergeCutHeight, corType = corType,
                          maxBlockSize = maxBlockSize, maxPOutliers = maxPOutliers,
                          numericLabels = numericLabels, pamRespectsDendro = pamRespectsDendro,
                          saveTOMs = TRUE,  saveTOMFileBase = dataset_name,
                          verbose = 0)
  
  moduleLabels = net$colors
  moduleColors = labels2colors(net$colors)
  MEs = net$MEs
  
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  moduleTraitCor = cor(MEs, traits, use = "p")
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  
  significant_module_colors <- sub(pattern = "ME", replacement = "", 
                                   x = rownames(moduleTraitPvalue)[stats::p.adjust(p = moduleTraitPvalue, method = "BH")
                                                                   < pval_cutoff])
  
  module_genes <- colnames(datExpr)[which(moduleColors %in% significant_module_colors)]
  
  info_table <- cbind(colnames(datExpr), moduleLabels, moduleColors)
  
  module_cor_and_p_value <- cbind(moduleTraitCor, moduleTraitPvalue, p.adjust(p = moduleTraitPvalue,
                                                                              method = "BH"))
  
  rownames(module_cor_and_p_value) <- sub(pattern = "ME", replacement = "", 
                                          x = rownames(module_cor_and_p_value))
  
  colnames(module_cor_and_p_value) <- c("Correlation", "p_value", "adjusted_p_value")
  
  new_wgcna_module <- wgcna_module_constructor(module_genes = module_genes,
                                               info_table = info_table,
                                               module_cor_and_p_value = module_cor_and_p_value,
                                               powerEstimate = powerEstimate,
                                               module_colors = significant_module_colors,
                                               input_class = class(MODifieR_input)[3],
                                               settings = settings)
  
  
  return(new_wgcna_module)
}

wgcna_get_trait_data <- function(MODifieR_input, group_of_interest){
  trait_data <- matrix(data = NA, ncol(MODifieR_input$annotated_exprs_matrix))
  if (group_of_interest == 1){
    trait_data[MODifieR_input$group_indici[[1]], ] <- 1
    trait_data[MODifieR_input$group_indici[[2]], ] <- 0
  }
  if (group_of_interest == 2){
    trait_data[MODifieR_input$group_indici[[1]], ] <- 0
    trait_data[MODifieR_input$group_indici[[2]], ] <- 1
  }
  rownames(trait_data) <- colnames(MODifieR_input$annotated_exprs_matrix)
  
  return (trait_data)
}
#' Generate a list of module colors with their respective genes. 
#'@param wgcna_module Module object that has been produced by \code{wgcna} function
#'@return Returns a \code{named list} of module genes where the names are the module colors. Includes 
#'non-significant colors 
#'
#'@seealso 
#'\code{\link{wgcna}}
#'
#'@author Dirk de Weerd
#'
#'
#'
wgcna_get_all_module_genes <- function(wgcna_module){
  module_colors <- rownames(wgcna_module$correlation_to_trait_table)
  
  module_genes_per_color <- sapply(X = module_colors, FUN = wgcna_get_module_genes, 
                                   info_table <- wgcna_module$info_table)
}

wgcna_get_module_genes <- function(module_color, info_table){
  info_table[info_table[,3] == module_color ,1 ]
}
#'@title Split WGCNA module in module containing only positive or negative correlation
#'@param wgcna_module Module object that has been produced by \code{wgcna} function
#'@param mode Character. "p" or "positive" for positive correlation, "n" or "negative" 
#'for negative correlation. 
#'@details
#' The functions returns a new \code{wgcna} module object that only contains positively or 
#' negatively correlated colors
#'
#' @seealso 
#' \code{\link{wgcna}}
#' 
#'@return \code{wgcna_module} object
#'
#'@author Dirk de Weerd
#'
#' @export
wgcna_get_module_genes_by_sign <- function(wgcna_module, mode){
  if (mode == "p" || mode =="positive"){
    module_colors <-  rownames(wgcna_module$correlation_to_trait_table)[
      wgcna_module$correlation_to_trait_table[, 1] > 0 
      & wgcna_module$correlation_to_trait_table[ ,3] < wgcna_module$settings$pval_cutoff]
  }else if (mode == "n" || mode == "negative"){
    module_colors <-  rownames(wgcna_module$correlation_to_trait_table)[
      wgcna_module$correlation_to_trait_table[, 1] < 0 
      & wgcna_module$correlation_to_trait_table[ ,3] < wgcna_module$settings$pval_cutoff]
  }else{
    warning('Incorrect mode argument provided. Use "postive" or "p" for postive correlation 
            #and "negative or "n" for negative correlation')
  }
  
  wgcna_module$module_genes <- wgcna_module$info_table[which(wgcna_module$info_table[ ,3] %in% module_colors), 1]
  
  return(wgcna_module)
}

#' wgcna_adjust_significance
#' @inheritParams clique_sum_permutation
#' @inheritParams wgcna_get_all_module_genes
#' @inheritParams wgcna
#' @param use_unadjusted Boolean value to signify if the adjusted (TRUE) or unadjusted (FALSE)
#' p value should be used to adjust significance
#' @details 
#' This function allows to adjust the significance cutoff for a \code{wgcna} module object
#' @return 
#'  \code{wgcna} module object 
#' @seealso 
#' 
#' \code{\link{wgcna}}
#' 
#' @author Dirk de Weerd
#' 
#' @export
wgcna_adjust_significance <- function(pval_cutoff, wgcna_module, use_unadjusted = F){
  col=3
  if (use_unadjusted){
    col=2
  }
  module_colors <- rownames(wgcna_module$correlation_to_trait_table)[
    wgcna_module$correlation_to_trait_table[ ,col] < pval_cutoff]
  
  wgcna_module$settings$pval_cutoff <- pval_cutoff
  wgcna_module$module_genes <- wgcna_module$info_table[which(wgcna_module$info_table[ ,3] %in% module_colors), 1]
  wgcna_module$module_colors <- module_colors
  
  return(wgcna_module)
}
#' Returns new WGCNA module objects by color
#' @inheritParams wgcna_get_all_module_genes
#' 
#' @details 
#' The  \code{wgcna} module object is split into a series of \code{wgcna} objects by color.
#' Every significant color in the module will be its own \code{wgcna} module object
#' 
#' @return 
#' 
#' A list of \code{wgcna} module objects
#' 
#' @seealso
#' 
#' \code{\link{wgcna}}
#' 
#' @author Dirk de Weerd
#' 
#' @export
wgcna_split_module_by_color <- function(wgcna_module){
  module_colors <- wgcna_module$module_colors
  module_genes <- lapply(X = module_colors, FUN = function(x, module){
    module$info_table[module$info_table[ ,3] == x, 1]}, module = wgcna_module)
  
  info_table <- wgcna_module$info_table
  correlation_to_trait_table <- wgcna_module$correlation_to_trait_table
  softthreshold_value <- wgcna_module$softthreshold_value
  settings <- wgcna_module$settings
  
  new_wgcna_modules <- list()
  
  for (color in 1:length(module_colors)){
    new_wgcna_modules[[color]] <- wgcna_module_constructor(module_genes = module_genes[[color]], 
                                                           info_table = info_table, 
                                                           module_cor_and_p_value = correlation_to_trait_table, 
                                                           powerEstimate = softthreshold_value, 
                                                           module_colors = module_colors[color], 
                                                           input_class = class(wgcna_module)[3],
                                                           settings = settings)
  }
  return(new_wgcna_modules)
}
#' Returns a wgcna module closest to \code{size}
#' 
#' @inheritParams wgcna_get_all_module_genes
#' @param size The desired size of the resulting module
#' 
#' @details 
#'  The function starts with the co-expression module (color) with the lowest 
#'  p-value and gradually adds more co-expression modules until the module will be 
#'  as close as possible to size \code{size} 
#'  
#'  @return 
#'  
#'  \code{wgcna} module object 
#'  
#' @seealso 
#' 
#' \code{\link{wgcna}}
#' 
#' @author Dirk de Weerd
#' 
#' @export
wgcna_set_module_size <- function(size, wgcna_module){
  #Get length of all modules
  module_lengths <- lengths(wgcna_get_all_module_genes(wgcna_module))
  #Order them according to adjusted p-value
  module_lengths <- module_lengths[order(wgcna_module$correlation_to_trait_table[, 3])]
  #Get the index of the last color to add to make the module size as close to "module_size"
  color_index <- which.min(abs(sapply(X = 1:length(module_lengths), FUN = function(x){
    size - sum(module_lengths[1:x])
  })))
  module_colors <- names(module_lengths)[1:color_index]
  new_wgcna_module  <-       wgcna_module_constructor(module_genes = wgcna_module$info_table
                                                     [which(wgcna_module$info_table[ ,3] %in% module_colors), 1], 
                                                     info_table = wgcna_module$info_table, 
                                                     module_cor_and_p_value = wgcna_module$correlation_to_trait_table, 
                                                     powerEstimate = wgcna_module$softthreshold_value, 
                                                     module_colors = module_colors, 
                                                     input_class = class(wgcna_module)[3],
                                                     settings = wgcna_module$settings)
  
  return(new_wgcna_module)
  
}


get_softthreshold_value <- function(datExpr, nSamples = nSamples){
  powers <- c(c(1:10), seq(from = 12, to=20, by=2))
  
  allowWGCNAThreads()
  
  sft <- WGCNA::pickSoftThreshold(t(datExpr), powerVector = powers, verbose = 0)
  
  powerEstimate <- sft$powerEstimate
  #Values taken from:
  #https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/faq.html
  if(is.na(powerEstimate)){
    if (nSamples < 20){
      powerEstimate <- 9
    }
    if (nSamples >= 20 && nSamples < 30){
      powerEstimate <- 8
    }
    if (nSamples >= 30 && nSamples < 40){
      powerEstimate <- 7
    }
    if (nSamples >= 40){
      powerEstimate <- 6
    }
  }
  return(powerEstimate)
}