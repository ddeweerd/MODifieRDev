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

validate_inputs <- function(settings){
  unit_interval_vars <- names(settings)[names(settings) %in% unit_intervals]
  for (unit in unit_interval_vars){
    unit <- settings[unit]
    validate_unitinterval(number = unlist(unit), varname = names(unit))
  }
  integer_vars <- names(settings)[names(settings) %in% integer_variables]
  for (integer_var in integer_vars){
    integer_var <- settings[integer_var]
    validate_integer(number = unlist(integer_var), varname = names(integer_var))
  }
  boolean_vars <- names(settings)[names(settings) %in% boolean_variables]
  for (boolean_var in boolean_vars){
    boolean_var <- settings[boolean_var]
    validate_boolean(boolean = unlist(boolean_var), varname = names(boolean_var))
  }
  
  numeric_vars <- names(settings)[names(settings) %in% numeric_variables]
  for (numeric_var in numeric_vars){
    numeric_var <- settings[numeric_var]
    validate_numeric(numeric_var = unlist(numeric_var), varname = names(numeric_var))
  }

}

validate_unitinterval <- function(number, varname){
  if (length(number) > 1){
    stop(varname, " is more than one value", call. = F)
  }
  if (!is.numeric(number)){
    stop(varname, " is not a number", call. = F)
  }
  if (!(number >= 0 & number <= 1)){
    stop(varname, " is not between 0 and 1", call. = F)
  }
}

validate_integer <- function(number, varname){
  if (length(number) > 1){
    stop(varname, " is more than one value", call. = F)
  }
  if (!is.numeric(number)){
    stop(varname, " is not a number", call. = F)
  }
  if (!isTRUE(number == as.integer(number))){
    stop(varname, " is not an integer", call. = F)
  }
}

validate_boolean <- function(boolean, varname){
  if (length(boolean) > 1){
    stop(varname, " is more than one value", call. = F)
  }
  if (!is.logical(boolean)){
    stop(varname, " is not logical", call. = F)
  }
}

validate_numeric <- function(numeric_var, varname){
  if (length(numeric_var) > 1){
    stop(varname, " is more than one value", call. = F)
  }
  if (!is.numeric(numeric_var)){
    stop(varname, " is not a numeric value", call. = F)
  }
}

validate_matrix <- function(matrix_var, varname){
  if (!is.matrix(matrix_var)){
    stop(varname, " is not a matrix", call. = F)
  }
}

check_diff_genes <- function(MODifieR_input, deg_cutoff = NULL){
  if (is.null(MODifieR_input$diff_genes)){
    stop("Differential expression missing (MODifieR_input$diff_genes)", call. = F)
  }
  if(!is.null(deg_cutoff)){
    if (sum(MODifieR_input$diff_genes$pvalue < deg_cutoff) == 0){
      stop("No differentially expressed genes below ", deg_cutoff)
    }
  }
}

check_expression_matrix <- function(MODifieR_input){
  if (is.null(MODifieR_input$annotated_exprs_matrix)){
    stop("Annotated expression matrix missing (MODifieR_input$annotated_exprs_matrix)", call. = F)
  }
}

validate_indici <- function(group1_indici, group2_indici, n_samples){
  n_indici <- length(c(group1_indici, group2_indici))
  if (n_indici != n_samples){
    stop("Different number of indici, ", n_indici, " given and ", n_samples, " samples", call. = F)
  }
  if (!all(unique(sort(as.integer(c(group1_indici, group2_indici)))) == 
           seq(from = 1, to = n_samples, by = 1))){
    stop("Group indici does not match all samples")
  }

}

