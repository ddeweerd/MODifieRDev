preprocess_object <- function(MODifieR_object){
  MODifieR_object <- prepare_export(MODifieR_object)
  MODifieR_object <- matrix_to_dataframe(MODifieR_object)
  
  return(MODifieR_object)
}

prepare_object <- function(MODifieR_object){
  remove_setting <- names(MODifieR_object$settings)[names(MODifieR_object$settings) %in% unexported_settings]
  
  for (setting in remove_setting) {
    MODifieR_object$settings[[which(names(MODifieR_object$settings) == setting)]] <- NULL
  }
  
  concatenate_rownames <- names(MODifieR_object)[names(MODifieR_object) %in% exported_rownames]
  
  for (rowname in concatenate_rownames) {
    rowname_index <- which(names(MODifieR_object) == rowname)
    data_table <- cbind(rownames(MODifieR_object[[rowname_index]]), MODifieR_object[[rowname_index]])
    if (class(data_table) == "data.frame"){
      colnames(data_table)[1] <- ""
    }
    MODifieR_object[[rowname_index]] <- data_table 
  }
  
  concatenate_names <- names(MODifieR_object)[names(MODifieR_object) %in% exported_names]
  
  for (name in concatenate_names) {
    name_index <- which(names(MODifieR_object) == name)
    data_table <- cbind(names(MODifieR_object[[name_index]]), MODifieR_object[[name_index]])
    if (class(data_table) == "data.frame"){
      colnames(data_table)[1] <- ""
    }
    MODifieR_object[[name_index]] <- data_table 
  }
  
  quote_args <- names(MODifieR_object$settings) %in% MODifieRDev:::unevaluated_args
  MODifieR_object$settings[quote_args] <- as.character(MODifieR_object$settings[quote_args])
  MODifieR_object$settings <- as.data.frame(MODifieR_object$settings, row.names = NULL)
  
  return(MODifieR_object)
}

matrix_to_dataframe <- function(MODifieR_object){
  matrix_indici <- which(sapply(MODifieR_object, class) == "matrix")
  for (matrix_index in matrix_indici){
    data_frame <- as.data.frame(MODifieR_object[[matrix_index]])
    if ( colnames(data_frame)[1] == "V1"){
      colnames(data_frame)[1] <- ""
    }
    MODifieR_object[[matrix_index]] <- data_frame
  }
  
  return(MODifieR_object)
}

prepare_export <- function(MODifieR_object){
  UseMethod("prepare_export", MODifieR_object)
}

prepare_export.MODifieR_input <- function(MODifieR_object){
  
  MODifieR_object$settings$group1_indici <- NULL
  MODifieR_object$settings$group2_indici <- NULL
  MODifieR_object$group_indici[[which.min(lengths(MODifieR_object$group_indici))]] <- 
    c(MODifieR_object$group_indici[[which.min(lengths(MODifieR_object$group_indici))]],
      rep(x = NA, abs(length(MODifieR_object$group_indici[[1]]) - length(MODifieR_object$group_indici[[2]]))))
  MODifieR_object$group_indici <- as.data.frame(MODifieR_object$group_indici)
  
  
  return(MODifieR_object)
}
prepare_export.Clique_Sum_exact <- function(MODifieR_object){
  MODifieR_object
}

prepare_export.Clique_Sum_permutation <- function(MODifieR_object){
  MODifieR_object
}

prepare_export.Correlation_clique <- function(MODifieR_object){
  MODifieR_object$frequency_table <- as.data.frame(MODifieR_object$frequency_table)
  colnames(MODifieR_object$frequency_table) <- c("Gene", "Frequency")
  
  return (MODifieR_object)
}

prepare_export.DIAMOnD <- function(MODifieR_object){
  MODifieR_object
}
prepare_export.DiffCoEx <- function(MODifieR_object){
  colnames(MODifieR_object$module_p_values) <- c("Color", "p value")
  color_frame <- MODifieR_object$color_vector[order(MODifieR_object$color_vector[,1]), ]
  colnames(color_frame) <- c("Color", "Gene")
  MODifieR_object$color_vector <- color_frame
  
  return(MODifieR_object)
}
prepare_export.Mcode <- function(MODifieR_object){
  module_lengths <- lengths(MODifieR_object$modules)
  module_names <- paste0("Module ", 1:length(module_lengths))
  module_frame <- get_module_frame(module_set = MODifieR_object$modules, module_lengths = module_lengths)
  rownames(module_frame) <- module_names
  MODifieR_object$modules <- t(module_frame)
  module_scores <- cbind(module_names, MODifieR_object$module_scores)
  colnames(module_scores) <- c("Module", "Score")
  MODifieR_object$module_scores <- module_scores
  
  return(MODifieR_object)
 
}

prepare_export.MODA<- function(MODifieR_object, xlsx_file){
  module_lengths1 <- lengths(MODifieR_object$group1_modules)
  module_lengths2 <- lengths(MODifieR_object$group2_modules)
  module_frame1 <- get_module_frame(module_set = MODifieR_object$group1_modules, module_lengths = module_lengths1)
  module_frame2 <- get_module_frame(module_set = MODifieR_object$group2_modules, module_lengths = module_lengths2)
  MODifieR_object$group1_modules <- t(module_frame1)
  MODifieR_object$group2_modules <- t(module_frame2)
  class(MODifieR_object) <- "list"
 
  return(MODifieR_object)
}

prepare_export.module_discoverer <- function(MODifieR_object, xlsx_file){
  MODifieR_object$graph <- NULL
  
  return(MODifieR_object)
}

prepare_export.WGCNA <- function(MODifieR_object, xlsx_file){
  colnames(MODifieR_object$info_table) <- c("Gene", "Label", "Color")
  
  return(MODifieR_object)
}

prepare_export.S2B <- function(MODifieR_object, xlsx_file){
  MODifieR_object
}

prepare_export.Module_set <- function(MODifieR_object, xlsx_file){
  module_gene_lengths <- lengths(MODifieR_object$module_gene_list)
  module_gene_frame <- t(get_module_frame(module_set = MODifieR_object$module_gene_list, 
                                          module_lengths = module_gene_lengths))
  colnames(module_gene_frame) <- names(MODifieR_object$module_gene_list)
  MODifieR_object$module_gene_list <- module_gene_frame
  gene_frequency <- as.data.frame(MODifieR_object$gene_frequency)
  colnames(gene_frequency) <- c("Gene", "Frequency")
  MODifieR_object$gene_frequency <- gene_frequency
  method_gene_lengths <- lengths(MODifieR_object$method_by_gene)
  method_frame <- t(get_module_frame(module_set = MODifieR_object$method_by_gene, 
                                     module_lengths = method_gene_lengths))
  colnames(method_frame)  <- names(MODifieR_object$method_by_gene)
  MODifieR_object$method_by_gene <- method_frame
 
  colnames(MODifieR_object$gene_by_method) <- c("Methods", "Frequency")
  
 
  return(MODifieR_object)
}
get_module_frame <- function(module_set, module_lengths){
  module_frame <- matrix(data = NA, nrow = length(module_lengths), ncol = max(module_lengths))
  for (i in 1:length(module_lengths)){
    module_frame[i, 1:module_lengths[i]] <- module_set[[i]]
  }
  return (module_frame)
}