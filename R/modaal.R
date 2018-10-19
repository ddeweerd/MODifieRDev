#' MODA
#'
#'A Coexpression based algorithm to identify disease module
#' @inheritParams clique_sum
#' @param datExpr1 A dataframe of the microarray matrix of the Control samples
#' @param datExpr2 A dataframe of the microarray matrix of the Control samples
#' @export
modaal <- function(MODifieR_input,
                   CuttingCriterion = "Density",
                   specificTheta = 0.1, conservedTheta = 0.1,
                   dataset_name = NULL, group_of_interest){
  # Retrieve settings
  default_args <- formals()
  user_args <- as.list(match.call(expand.dots = T)[-1])
  settings <- c(user_args, default_args[!names(default_args) %in% names(user_args)])
  
  if (!is.null(dataset_name)){
    settings$MODifieR_input <- dataset_name
  }
  #Get relevant input data from input object 
  datExpr1 <- t(MODifieR_input$annotated_exprs_matrix[,MODifieR_input$group_indici[[1]]])
  datExpr2 <- t(MODifieR_input$annotated_exprs_matrix[,MODifieR_input$group_indici[[2]]])
  
  indicator1 = names(MODifieR_input$group_indici)[1]     # indicator for data profile 1
  indicator2 = names(MODifieR_input$group_indici)[2]  # indicator for data profile 2
  
  modules1 <- WMPH(datExpr = datExpr1, indicatename = indicator1, cutmethod = CuttingCriterion)
  
  modules2 <- WMPH(datExpr = datExpr2, indicatename = indicator2, cutmethod = CuttingCriterion)
  
  if (group_of_interest == 1){
    module_score_array <- compare_modules(module_list1 = modules1, module_list2 = modules2)
    
    module_genes <- unique(unlist(modules1[moda_extract_modules_index_specific(
      module_score_array = module_score_array, specificTheta = specificTheta)]))
  }
  if (group_of_interest == 2){
    module_score_array <- compare_modules(module_list1 = modules2, module_list2 = modules1)
    
    module_genes <- unique(unlist(modules2[moda_extract_modules_index_specific(
      module_score_array = module_score_array, specificTheta = specificTheta)]))
  }
  
  new_moda_module <- list("module_genes" =  module_genes,
                          "group1_modules" = modules1,
                          "group2_modules" = modules2,
                          "module_score_array" = module_score_array,
                          "settings" = settings)
  
  class(new_moda_module) <- c("MODifieR_module", "MODA")
  
  return (new_moda_module)
}
