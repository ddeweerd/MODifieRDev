validate_new_input_objects <- function(...){
  
  arg_list  <- list(...)
  if (!is.null(arg_list$diff_genes)){
    arg_list$diff_genes <- validate_diff_genes(diff_genes = arg_list$diff_genes)
  }
  if (!is.null(arg_list$annotated_exprs_matrix)){
    
    validated_expression <- validate_exprs_matrix(annotated_exprs_matrix = arg_list$annotated_exprs_matrix, 
                                                    group1_indici = arg_list$group1_indici, 
                                                    group2_indici = arg_list$group2_indici, 
                                                    group1_label = arg_list$group1_label, 
                                                    group2_label = arg_list$group2_label)
    arg_list$annotated_exprs_matrix <- validated_expression$annotated_exprs_matrix
    arg_list$group1_indici <- validated_expression$group1_indici
    arg_list$group2_indici <- validated_expression$group2_indici
    arg_list$group1_label <- validated_expression$group1_label
    arg_list$group2_label <- validated_expression$group2_label
    
  }
  return(list("diff_genes" = arg_list$diff_genes,
              "annotated_exprs_matrix" = arg_list$annotated_exprs_matrix,
              "group1_indici" = arg_list$group1_indici,
              "group2_indici" = arg_list$group2_indici,
              "group1_label" =  arg_list$group1_label,
              "group2_label" = arg_list$group2_label))
}

validate_diff_genes <- function(diff_genes){
  if (class(diff_genes) == "data.frame"){
    if (ncol(diff_genes) != 2){
      stop("diff_genes has a different number of columns: ", ncol(diff_genes), call. = F)
    }
  }else{
    stop("diff_genes is not a data.frame", call. = F)
  }
  diff_genes[,1] <- as.character(diff_genes[,1])
  diff_genes[,2] <- as.numeric(diff_genes[,2])
  
  colnames(diff_genes) <- c("gene", "pvalue")
  
  return(diff_genes)
}

validate_exprs_matrix <- function(annotated_exprs_matrix, group1_indici,
                                  group2_indici, group1_label, group2_label){
  evaluated_args <- c(as.list(environment()))
  
  null_args <- names(evaluated_args[sapply(evaluated_args, is.null)])
  if (length(null_args > 0)){
    stop("Missing values: ", sapply(null_args, function(x){c(x, " ")}), call. = F)
  }
  if (class(annotated_exprs_matrix) != "matrix"){
    stop("annotated_exprs_matrix is not a matrix", call. = F)
  }
  return(list("annotated_exprs_matrix" = annotated_exprs_matrix,
              "group1_indici" = group1_indici,
              "group2_indici" = group2_indici,
              "group1_label" =  group1_label,
              "group2_label" = group2_label))
  
}

validate_ppi <- function(ppi_network){
  if (class(ppi_network) == "data.frame"){
    if (!ncol(ppi_network) %in% 2:3 ){
      stop("ppi_network has a different number of columns: ", ncol(ppi_network), call. = F)
    }
  }else{
    stop("ppi_network is not a data.frame", call. = F)
  }
  ppi_network[,1] <- as.character(ppi_network[,1])
  ppi_network[,2] <- as.character(ppi_network[,2])
  if (ncol(ppi_network) == 3){
    ppi_network[,1] <- as.numeric(ppi_network[,3])
  }
  return(ppi_network)
}
