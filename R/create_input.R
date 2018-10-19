#' Creates an input object for downstream analysis 
#' 
#' The MODifieR input object can be used in downstream analysis for the disease module inference methods included
#' in this package.
#'@param expression_matrix Normalized expression matrix where the samples are columns and probes are rows
#'@param annotation_table A dataframe providing annotation for the probes. The dataframe should have 3 columns: 
#' \itemize{
#' \item {PROBEID}: The probe id as it is in the expression matrix
#' \item {SYMBOL}: The gene symbol (if available) associated with the probe
#' \item{ENTREZID}: The entrez id (if available) associated with the probe
#'}
#'@param group_indici vectors containing indici for different groups (Column numbers)
#'@param group_labels Labels for each group, for example "patient" and "control"
#'@param expression boolean, calculate expression values?
#'@param differential_expression boolean, calculate differentially expressed data?
#'@param correlation_clique boolean, calculate correlation matrix?
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
#' \item{diff_genes}{A 2 two column data.frame where the first column are genes and the second column unadjusted p-values}
#' \item{limma_probe_table}{A data.frame from \code{limma topTable} with added gene annotation}
#' \item{annotated_exprs_matrix}{A matrix where the rows are genes and the columns samples. Probes have been collapsed
#' into genes using \code{collapse_method}}
#' \item{expression_matrix}{A matrix, the original input expression matrix}
#' \item{annotation_table}{A data.frame, the original annotation table used to annotate the probes}
#' \item{group_indici}{A named list containg 2 numeric vectors. The names are the group labels and the values 
#' are the group indici}
#' @export
create_input <- function (expression_matrix, probe_map, group1_indici, group2_indici, group1_label, group2_label,
                          expression = T,  differential_expression= T, method = "MaxMean"){
  #Initialize outputs
  diff_genes <- NULL
  collapsed_exprs_mat <- NULL
  limma_probe_table <- NULL
  group_indici <- list("group_1_indici" = group1_indici,
                       "group_2_indici" = group2_indici)
  names(group_indici) <- c(group1_label, group2_label)
  
  #Making sure column names for the annotation dataframe are right...
  colnames(probe_map) <- c("PROBEID", "SYMBOL", "ENTREZID")
  
  #Should expression data be generated? 
  if (expression == T){
    collapsed_data <-WGCNA::collapseRows(datET = expression_matrix,
                                         rowGroup =  probe_map$ENTREZID,
                                         rowID =  probe_map$PROBEID, method = collapse_method)
    
    collapsed_exprs_mat <- collapsed_data$datETcollapsed
  }
  #Same for expression data
  if (diff_data == T){
    group_factor <- create_group_factor(samples = colnames(expression_matrix),
                                        group1_indici = group1_indici,
                                        group2_indici = group2_indici)
    
    limma_probe_table <- differential_expression(group_factor = group_factor,
                                               expression_matrix = expression_matrix,
                                               probe_table = probe_map)
    
    diff_genes <- plyr::ddply(.data = , limma_probe_table,
                             .variables = "ENTREZID", .fun = plyr::summarise, pvalue = min(P.Value))
    
    diff_genes <- na.omit(diff_genes)
  }
  modifier_input <- list("diff_genes" = diff_genes,
                         "limma_probe_table" = limma_probe_table,
                         "annotated_exprs_matrix" = collapsed_exprs_mat,
                         "expression_matrix" = expression_matrix,
                         "annotation_table" = probe_map,
                         "group_indici" = group_indici)
  class(modifier_input) <- c("MODifieR_input", "Expression")
  return (modifier_input)
}
#Calculate differentially expressed genes
differential_expression <- function(group_factor, expression_matrix, probe_table){
  design <- model.matrix(~group_factor)
  fit <- limma::lmFit(expression_matrix, design)
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