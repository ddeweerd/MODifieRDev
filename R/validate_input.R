validate_input <- function(input_object){
  
  if(is.null(input_object$diff_genes)){
    warning("No differentially expression data\n")
  }
  if(is.null(input_object$annotated_exprs_matrix)){
    warning("No collapsed expression data\n")
  }
  if(is.null(input_object$expression_matrix)){
    warning("No expression matrix\n")
  }
  if(is.null(input_object$diamond_genes)){
    warning("No diamond genes\n")
  }
  if(is.null(input_object$annotation_table)){
    warning("No probe annotation table\n")
  }
  if(is.null(input_object$group_indici)){
    warning("No group indici\n")
  }
  if(length(unlist(input_object$group_indici)) != ncol(input_object$expression_matrix)){
    warning("Group indici do not match the number of columns in the expression matrix\n")
  }
  if(length(unlist(input_object$group_indici)) != ncol(input_object$annotated_exprs_matrix)){
    warning("Group indici do not match the number of columns in the collapsed expression matrix\n")
  }
  if(ncol(input_object$expression_matrix) != ncol(input_object$expression_matrix)){
    warning("Number of columns differ between collapsed expression matrix and expression matrix\n")
  }
  n_na_expr_mat <- sum(is.na(input_object$expression_matrix))
  n_na_cola_mat <- sum(is.na(input_object$annotated_exprs_matrix))
  
  if (n_na_expr_mat > 0){
    warning(paste("Expression matrix contains ", n_na_expr_mat, " NA values\n"))
  }
  if (n_na_cola_mat > 0){
    warning(paste("Collapsed expression matrix contains ", n_na_cola_mat, " NA values\n"))
  }
}

