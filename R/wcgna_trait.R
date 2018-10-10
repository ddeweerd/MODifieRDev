#Module object constructor
wgcna_module_constructor <- function(module_genes, probe_info_table, 
                                     module_cor_and_p_value, powerEstimate,
                                     module_colors, settings){
  
  new_wgcna_module <- list("module_genes" =  module_genes,
                           "probe_info_table" = probe_info_table,
                           "correlation_to_trait_table" = module_cor_and_p_value,
                           "softthreshold_value" = powerEstimate,
                           "module_colors" = module_colors,
                           "settings" = settings)
  
  class(new_wgcna_module) <- c("MODifieR_module", "WGCNA")
  
  return(new_wgcna_module)
}

#' An implementation of WGCNA to correlate coexpression modules to disease
#'
#' @param MODifieR_input An input object created by \code{\link[WGCNA]{blockwiseModules}} function
#' Paap: \code{\link[WGCNA]{blockwiseModules}}
#' @seealso \code{\link[MASS]{abbey}}
#' @export
wgcna_trait <- function(MODifieR_input,  minModuleSize = 30, deepSplit = 2, pamRespectsDendro = F,
                        mergeCutHeight = 0.1, numericLabels = T, min_KME = 0, pval_cutoff = 0.05,
                        save_toms = T, dataset_name = deparse(substitute(MODifieR_input))){
  # Retrieve settings
  default_args <- formals()
  user_args <- as.list(match.call(expand.dots = T)[-1])
  settings <- c(user_args, default_args[!names(default_args) %in% names(user_args)])
  settings$dataset_name <- NULL
  
  if (!is.null(dataset_name)){
    settings$MODifieR_input <- dataset_name
  }
  
  allowWGCNAThreads()
  
  traits <- wgcna_get_trait_data(MODifieR_input = MODifieR_input)
  
  datExpr <- t(MODifieR_input$annotated_exprs_matrix)
  
  nGenes <- ncol(datExpr)
  nSamples <- nrow(datExpr)
  
  powers <- c(c(1:10), seq(from = 12, to=20, by=2))
  
  sft <- pickSoftThreshold(t(datExpr), powerVector = powers, verbose = 0)
  
  powerEstimate <- sft$powerEstimate
  #Values taken from:
  #https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/faq.html
  if(is.na(powerEstimate)){
    if (nSamples < 20){
      powerEstimate <- 6
    }
    if (nSamples >= 20 && nSamples < 30){
      powerEstimate <- 7
    }
    if (nSamples >= 30 && nSamples < 40){
      powerEstimate <- 8
    }
    if (nSamples >= 40){
      powerEstimate <- 9
    }
  }
  net <- blockwiseModules(datExpr, power = powerEstimate,
                          TOMType = "unsigned", minModuleSize = minModuleSize,
                          mergeCutHeight = mergeCutHeight,
                          numericLabels = numericLabels, pamRespectsDendro = pamRespectsDendro,
                          saveTOMs = TRUE,  saveTOMFileBase = dataset_name,
                          verbose = 5)
  
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
                                   x = rownames(moduleTraitPvalue)[p.adjust(p = moduleTraitPvalue, method = "BH")
                                                                   < pval_cutoff])
  
  module_genes <- colnames(datExpr)[which(moduleColors %in% significant_module_colors)]
  
  probe_info_table <- cbind(colnames(datExpr), moduleLabels, moduleColors)
  
  module_cor_and_p_value <- cbind(moduleTraitCor, moduleTraitPvalue, p.adjust(p = moduleTraitPvalue,
                                                                              method = "BH"))
  
  rownames(module_cor_and_p_value) <- sub(pattern = "ME", replacement = "", 
                                          x = rownames(module_cor_and_p_value))
  
  colnames(module_cor_and_p_value) <- c("Correlation", "p_value", "adjusted_p_value")
  
  new_wgcna_module <- wgcna_module_constructor(module_genes = module_genes,
                                               probe_info_table = probe_info_table,
                                               module_cor_and_p_value = module_cor_and_p_value,
                                               powerEstimate = powerEstimate,
                                               module_colors = significant_module_colors,
                                               settings = settings)
  
  return(new_wgcna_module)
}

wgcna_get_trait_data <- function(MODifieR_input){
  trait_data <- matrix(data = NA, ncol(MODifieR_input$annotated_exprs_matrix))
  trait_data[MODifieR_input$group_indici[[1]], ] <- 0
  trait_data[MODifieR_input$group_indici[[2]], ] <- 1
  rownames(trait_data) <- colnames(MODifieR_input$annotated_exprs_matrix)
  
  return (trait_data)
}
#' Generate a list of module colors with their respective genes. 
#'@param wgcna_module Module object that has been produced by \code{wgcna_trait} function
#'@return Returns a \code{named list} of module genes where the names are the module colors. Includes 
#'non-significant colors 
#'@export
wgcna_get_all_module_genes <- function(wgcna_module){
  module_colors <- rownames(wgcna_module$correlation_to_trait_table)
  
  module_genes_per_color <- sapply(X = module_colors, FUN = wgcna_get_module_genes, 
                                   probe_info_table <- wgcna_module$probe_info_table)
}

wgcna_get_module_genes <- function(module_color, probe_info_table){
  probe_info_table[probe_info_table[,3] == module_color ,1 ]
}
#'@title Split WGCNA module in module containing only positive or negative correlation
#'@param wgcna_module Module object that has been produced by \code{wgcna_trait} function
#'@param mode Character. "p" or "positive" for positive correlation, "n" or "negative" 
#'for negative correlation
#'
#'Only significant colors are used
#'@return Returns a \code{wgcna_module} object with only postively/negatively associated module colors to the trait
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
  
  wgcna_module$module_genes <- wgcna_module$probe_info_table[which(wgcna_module$probe_info_table[ ,3] %in% module_colors), 1]
  
  return(wgcna_module)
}
#' @export
wgcna_adjust_significance <- function(p_value, wgcna_module, use_unadjusted = F){
  col=3
  if (use_unadjusted){
    col=2
  }
  module_colors <- rownames(wgcna_module$correlation_to_trait_table)[
  wgcna_module$correlation_to_trait_table[ ,col] < p_value]
  
  wgcna_module$settings$pval_cutoff <- p_value
  wgcna_module$module_genes <- wgcna_module$probe_info_table[which(wgcna_module$probe_info_table[ ,3] %in% module_colors), 1]
  
  return(wgcna_module)
}
#Returns new module objects by color
#' @export
wgcna_split_module_by_color <- function(wgcna_module){
  module_colors <- wgcna_module$module_colors
  module_genes <- lapply(X = module_colors, FUN = function(x, module){
  module$probe_info_table[module$probe_info_table[ ,3] == x, 1]}, module = wgcna_module)
  
  probe_info_table <- wgcna_module$probe_info_table
  correlation_to_trait_table <- wgcna_module$correlation_to_trait_table
  softthreshold_value <- wgcna_module$softthreshold_value
  settings <- wgcna_module$settings
  
  new_wgcna_modules <- list()
  
  for (color in 1:length(module_colors)){
    new_wgcna_modules[[color]] <- wgcna_module_constructor(module_genes = module_genes[[color]], 
                                                           probe_info_table = probe_info_table, 
                                                           module_cor_and_p_value = correlation_to_trait_table, 
                                                           powerEstimate = softthreshold_value, 
                                                           module_colors = module_colors[color], 
                                                           settings = settings)
  }
  return(new_wgcna_modules)
}
#' Returns a wgcna module close to \code{parameter}
#' @export
wgcna_set_module_size <- function(size, wgcna_module){
  counter <- 0
  module_colors <- NULL
  while (size > 0){
    counter <- counter + 1
    current_color <- lengths(wgcna_get_all_module_genes(wgcna_module))[(order(wgcna_module$correlation_to_trait_table[,2]))][counter]
    size <- size - current_color
    module_colors[counter] <- names(current_color)
    
  }
  new_wgcna_module <-      wgcna_module_constructor(module_genes = wgcna_module$probe_info_table
                                                    [which(wgcna_module$probe_info_table[ ,3] %in% module_colors), 1], 
                                                    probe_info_table = wgcna_module$probe_info_table, 
                                                    module_cor_and_p_value = wgcna_module$correlation_to_trait_table, 
                                                    powerEstimate = wgcna_module$softthreshold_value, 
                                                    module_colors = module_colors, 
                                                    settings = wgcna_module$settings)
  
  return(new_wgcna_module)
  
}
