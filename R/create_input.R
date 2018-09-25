#' Creates an object of class "MODifieR input" 
#' 
#' The MODifieR input object can be used
#' in downstream analysis using the disease module inference methods included
#' in the MODifieR package.
#'@param expression_matrix Normalized expression matrix where the samples are columns and probes are rows
#'@param probe_map A dataframe providing annotation for the probes. The dataframe should have 3 columns: 
#' \itemize{
#' \item {PROBEID}: The probe id as it is in the expression matrix
#' \item {SYMBOL}: The gene symbol (if available) associated with the probe
#' \item{ENTREZID}: The entrez id (if available) associated with the probe
#'}
#'@param group_indici vectors containing indici for different groups (Column numbers)
#'@param group_labels Labels for each group, for example "patient" and "control"
#'@param expression boolean, calculate expression values?
#'@param diff_data boolean, calculate differentially expressed data?
#'@param correlation_clique boolean, calculate correlation matrix?
#' 
#' @export
create_input <- function (expression_matrix, probe_map, group1_indici, group2_indici, group1_label, group2_label,
                             expression = F, diff_data = F, collapse_method = "MaxMean"){
  #Initialize outputs
  diff_genes <- NULL
  collapsed_exprs_mat <- NULL
  diamond_genes <- NULL
  group_indici <- list("group_1_indici" = group1_indici,
                       "group_2_indici" = group2_indici)
  names(group_indici) <- c(group1_label, group2_label)

  #Making sure column names for the annotation dataframe are right...
  colnames(probe_map) <- c("PROBEID", "SYMBOL", "ENTREZID")

  #Should expression data be generated? Note that expression data is precursor
  #for correlation method data
  if (expression == T){
    collapsed_data <-WGCNA::collapseRows(datET = expression_matrix,
                                         rowGroup =  probe_map$ENTREZID,
                                         rowID =  probe_map$PROBEID, method = collapse_method)
    
    collapsed_exprs_mat <- collapsed_data$datETcollapsed
    
  }
  #Same for expression data, which is also a precursor to correlation data.
  if (diff_data == T){
    group_factor <- create_group_factor(samples = colnames(expression_matrix),
                                         group1_indici = group1_indici,
                                         group2_indici = group2_indici)

    diff_genes_data <- differential_expression(group_factor = group_factor,
                                                expression_matrix = expression_matrix,
                                                probe_table = probe_map)
    if (diff_data == T){
      diff_genes <- data.frame(gene = diff_genes_data$ENTREZID ,
                               pvalue = diff_genes_data$P.Value,
                               stringsAsFactors = FALSE)
      diff_genes <- na.omit(diff_genes)
      diff_genes <- summarise_pval(pvalues =  diff_genes)
      diamond_genes <- diff_genes_data[diff_genes_data$P.Value < 0.05, ]
      diamond_genes <- unique(na.omit(diamond_genes$ENTREZID)) 
      }
  }
  modifier_input <- list("diff_genes" = diff_genes,
                         "annotated_exprs_matrix" = collapsed_exprs_mat,
                         "expression_matrix" = expression_matrix,
                         "diamond_genes" = diamond_genes,
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
#Helper function for correlation clique
summarise_pval <- function(pvalues){
  colnames(pvalues) <- c("gene", "pvalue")
  plyr::ddply(.data = pvalues, .variables = "gene", .fun = plyr::summarise, pvalue = min(pvalue))
}
#Collapse probes that target same gene and get mean expression value.
collapse_probes <- function(entrez, exprs_mat, probe_entrez_mat){
  sapply(X = 1:ncol(exprs_mat),
         FUN = function(x){mean(exprs_mat[probe_entrez_mat[which(probe_entrez_mat$ENTREZID == entrez),1], x])})
}
