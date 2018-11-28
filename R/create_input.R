#' Creates an input object for downstream analysis 
#' 
#' The MODifieR input object can be used in downstream analysis for the disease module inference methods included
#' in this package.
#'@param expression_matrix Normalized expression matrix where the samples are columns and probes are rows
#'@param annotation_table A dataframe providing annotation for the probes. The dataframe should have 3 columns: 
#' \itemize{
#' \item {PROBEID}: The probe id as it is in the expression matrix
#' \item{IDENTIFIER}: The entrez id (if available) associated with the probe
#'}
#'@param group1_indici vector containing indici for samples belonging to group 1 (Column numbers)
#'@param group2_indici vector containing indici for samples belonging to group 2 (Column numbers)
#'@param group1_label Label for each group 1, for example "patient" or "control"
#'@param group2_label Label for each group 2, for example "patient" or "control"
#'@param expression boolean, calculate expression values?
#'@param differential_expression boolean, calculate differentially expressed data?
#' @inheritParams WGCNA::collapseRows
#' @details 
#' The function creates an input object to be used in all disease module inference methods. Differentially
#' expressed genes are calculated using linear models from the \code{limma} package. Probes are collapsed into genes
#' using collapseRows from \code{WGCNA}
#' @seealso 
#' \code{\link[WGCNA]{collapseRows}}
#' \code{\link[limma]{lmFit}}
#' \code{\link[limma]{eBayes}}
#' \code{\link[limma]{topTable}}
#' @return
#' The function returns an object of class "MODifieR_input". The object is a named list containing the
#' following components:
#' \item{diff_genes}{A 2 two column data.frame where the first column are genes and the second column unadjusted p-values
#' obtained by differential expression analysis}
#' \item{limma_probe_table}{A data.frame from \code{limma topTable} with added gene annotation}
#' \item{annotated_exprs_matrix}{A matrix where the rows are genes and the columns samples. Probes have been collapsed
#' into genes using \code{collapse_method}}
#' \item{expression_matrix}{A matrix, the original input expression matrix}
#' \item{annotation_table}{A data.frame, the original annotation table used to annotate the probes}
#' \item{group_indici}{A named list containing 2 numeric vectors. The names are the group labels and the values 
#' are the group indici}
#' @export
create_input <- function (expression_matrix, annotation_table, group1_indici, group2_indici, group1_label, group2_label,
                          expression = T,  differential_expression= T, method = "MaxMean"){
  #Initialize outputs
  diff_genes <- NULL
  collapsed_exprs_mat <- NULL
  limma_probe_table <- NULL

  #Making sure column names for the annotation dataframe are right...
  colnames(annotation_table) <- c("PROBEID", "IDENTIFIER")
  
  #Should expression data be generated? 
  if (expression == T){
    
    freq_probes <- table(annotation_table$PROBEID)
    
    probe_names <- names(freq_probes[freq_probes == 1])
    
    annotation_table_expression <- annotation_table[annotation_table$PROBEID  %in% probe_names,]
    
    annotation_table_expression <- stats::na.omit(annotation_table_expression)
    
    collapsed_data <-WGCNA::collapseRows(datET = expression_matrix,
                                         rowGroup =  annotation_table_expression$IDENTIFIER,
                                         rowID =  annotation_table_expression$PROBEID, method = method)
    
    collapsed_exprs_mat <- collapsed_data$datETcollapsed
  }
  #Same for expression data
  if (differential_expression == T){
    P.Value = NULL
    group_factor <- create_group_factor(samples = colnames(expression_matrix),
                                        group1_indici = group1_indici,
                                        group2_indici = group2_indici)
    
    limma_probe_table <- differential_expression(group_factor = group_factor,
                                               expression_matrix = expression_matrix,
                                               probe_table = annotation_table)
    
    diff_genes <- plyr::ddply(.data = , limma_probe_table,
                             .variables = "IDENTIFIER", .fun = plyr::summarise, pvalue = min(P.Value))
    
    diff_genes <- stats::na.omit(diff_genes)
    colnames(diff_genes) <- c("gene", "pvalue")
  }
  modifier_input <- create_custom_input_object(diff_genes = diff_genes, 
                                           limma_probe_table = limma_probe_table,
                                           annotated_exprs_matrix = collapsed_exprs_mat,
                                           expression_matrix = expression_matrix,
                                           group1_indici = group1_indici,
                                           group2_indici = group2_indici,
                                           group1_label = group1_label,
                                           group2_label = group2_label)
 }
#Calculate differentially expressed genes
differential_expression <- function(group_factor, expression_matrix, probe_table){
  design <- stats::model.matrix(~group_factor)
  fit <- invisible(limma::lmFit(expression_matrix, design))
  fit2 <- limma::eBayes(fit)
  diff_genes <- limma::topTable(fit = fit2,  number = Inf, adjust.method = "BH")
  annotated_probes <- lapply(X = rownames(diff_genes), FUN = annotate_probe,
                             diff_genes = diff_genes, probe_table = probe_table)
  annotated_probes <- do.call("rbind", annotated_probes)
  rownames(annotated_probes) <- NULL
  return (annotated_probes)
}
#Helper function for differential expression, adds annotation to probes
annotate_probe <- function(probe_id, diff_genes, probe_table){
  cbind(diff_genes[probe_id, ], probe_table[probe_table$PROBEID == probe_id, ])
}
#Helper function for differential expression, creates group factor
create_group_factor <- function(samples, group1_indici, group2_indici){
  pre_factor <- samples
  pre_factor <- replace(x = pre_factor, list = group1_indici, values = "group1")
  pre_factor <- replace(x = pre_factor, list = group2_indici, values = "group2")
  group_factor<- as.factor(pre_factor)
  return(group_factor)
}
#' Create a generic input object
#' @inheritParams create_input
#' @param diff_genes A 2 two column data.frame where the first column are genes and the second column are p-values
#' @param limma_probe_table A data.frame from \code{limma topTable} with added gene annotation
#' @param annotated_exprs_matrix A matrix where the rows are genes and the columns samples.
#' 
#' @details 
#' This function allows the creation of a generic input object with the same class as objects created by
#' \code{\link{create_input}}. This can be useful in the cases where you already have differentially 
#' expressed genes, an annotated expression matrix or both and want to wrap that into an input object
#' to use in downstream analysis.
#' @seealso 
#' \code{\link{create_input}}
#' @return
#' The function returns an object of class "MODifieR_input". The object is a named list containing the
#' following components:
#' \item{diff_genes}{A 2 two column data.frame where the first column are genes and the second column unadjusted p-values}
#' \item{limma_probe_table}{A data.frame from \code{limma topTable} with added gene annotation}
#' \item{annotated_exprs_matrix}{A matrix where the rows are genes and the columns samples. Probes have been collapsed
#' into genes using \code{collapse_method}}
#' \item{expression_matrix}{A matrix, the original input expression matrix}
#' \item{annotation_table}{A data.frame, the original annotation table used to annotate the probes}
#' \item{group_indici}{A named list containing 2 numeric vectors. The names are the group labels and the values 
#' are the group indici}
create_custom_input_object <- function(diff_genes = NULL, limma_probe_table = NULL,
                                   annotated_exprs_matrix = NULL, expression_matrix = NULL, 
                                   annotation_table = NULL, group1_indici = NULL,
                                   group2_indici = NULL, group1_label = NULL, 
                                   group2_label = NULL){
  
  group_indici <- list("group_1_indici" = group1_indici,
                       "group_2_indici" = group2_indici)
  names(group_indici) <- c(group1_label, group2_label)
  
  modifier_input <- list("diff_genes" = diff_genes,
                         "limma_probe_table" = limma_probe_table,
                         "annotated_exprs_matrix" = annotated_exprs_matrix,
                         "expression_matrix" = expression_matrix,
                         "annotation_table" = annotation_table,
                         "group_indici" = group_indici)
  class(modifier_input) <- c("MODifieR_input", "Expression")
  
  return (modifier_input)
  
  
}
