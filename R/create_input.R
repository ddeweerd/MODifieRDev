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
#'@param filter_expression boolean, remove 50 percent of the genes with lowest variance?
#'@param use_adjusted boolean, use adjusted p value for differential expression analysis?
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
#' 
#' @author Dirk de Weerd
#' @export
create_input_microarray <- function (expression_matrix, annotation_table, group1_indici, group2_indici, group1_label, group2_label,
                                     expression = T,  differential_expression= T, method = "MaxMean", filter_expression = T,
                                     use_adjusted = T){
  
  # Retrieve settings
  evaluated_args <- c(as.list(environment()))
  settings <- as.list(stackoverflow::match.call.defaults()[-1])
  replace_args <- names(settings)[!names(settings) %in% unevaluated_args]
  for (argument in replace_args) {
    settings[[which(names(settings) == argument)]] <- evaluated_args[[which(names(evaluated_args) == 
                                                                              argument)]]
  }
  #Initialize outputs
  diff_genes <- NULL
  collapsed_exprs_mat <- NULL
  limma_probe_table <- NULL
  if (filter_expression){
    exp_matrix_var <- apply(expression_matrix, 1, function(x) stats::var(x, na.rm = TRUE))
    expression_matrix <- expression_matrix[exp_matrix_var >= stats::quantile(exp_matrix_var, c(.50)), ]
  }
  
  #Making sure column names for the annotation dataframe are right...
  colnames(annotation_table) <- c("PROBEID", "IDENTIFIER")
  #And of the right class
  annotation_table$IDENTIFIER <- as.character(annotation_table$IDENTIFIER)
  #Should expression data be generated? 
  if (expression == T){
    collapsed_exprs_mat <- collapse_probes(expression_matrix = expression_matrix, 
                                           annotation_table = annotation_table, 
                                           method = method)
  }
  #Same for expression data
  if (differential_expression == T){
    P.Value = NULL
    adj.P.Val = NULL
    
    group_factor <- create_group_factor(samples = colnames(expression_matrix),
                                        group1_indici = group1_indici,
                                        group2_indici = group2_indici)
    
    limma_probe_table <- differential_expression(group_factor = group_factor,
                                                 expression_matrix = expression_matrix,
                                                 probe_table = annotation_table)
    
    diff_genes <- summarize_probes(limma_probe_table = limma_probe_table, use_adjusted = use_adjusted)
    
    diff_genes <- stats::na.omit(diff_genes)
    
    
  }
  modifier_input <- create_custom_microarray_input_object(diff_genes = diff_genes, 
                                                          limma_probe_table = limma_probe_table,
                                                          annotated_exprs_matrix = collapsed_exprs_mat,
                                                          expression_matrix = expression_matrix,
                                                          annotation_table = annotation_table,
                                                          group1_indici = group1_indici,
                                                          group2_indici = group2_indici,
                                                          group1_label = group1_label,
                                                          group2_label = group2_label,
                                                          settings)
  
  return (modifier_input)
}


#' Creates an RNA-seq input object for downstream analysis 
#' 
#' The MODifieR input objects can be used in downstream analysis for the disease module inference methods included
#' in this package.
#' @param count_matrix Matrix containing raw RNA-seq counts
#' @inheritParams create_input_microarray
#' @param normalize_quantiles boolean, Normalize quantiles for WGCNA-based methods?
#' @details 
#' The function creates an input object to be used in all disease module inference methods. Differentially
#' expressed genes are calculated using generalized linear models from the \code{edgeR} package. For WGCNA-based
#' methods raw counts are normalized using the \code{varianceStabilizingTransformation} from the \code{DESeq2} package. 
#' Optionally, Quantile normalization using the \code{normalize.quantiles} function from the \code{preprocessCore}
#' can be applied.
#' @seealso 
#' \code{\link[edgeR]{glmQLFit}}
#' \code{\link[DESeq2]{varianceStabilizingTransformation}}
#' @return
#' The function returns an object of class "MODifieR_input". The object is a named list containing the
#' following components:
#' \item{diff_genes}{A 2 two column data.frame where the first column are genes and the second column unadjusted p-values
#' obtained by differential expression analysis}
#' \item{edgeR_deg_table}{A data.frame from \code{edgeR glmQLFit}}
#' \item{annotated_exprs_matrix}{A matrix where the rows are genes and the columns samples. Raw counts have been normalized
#' with \code{varianceStabilizingTransformation} and optionally quantile normalized}
#' \item{count_matrix}{A matrix, the original input count matrix}
#' \item{group_indici}{A named list containing 2 numeric vectors. The names are the group labels and the values 
#' are the group indici}
#' @export
create_input_rnaseq <- function(count_matrix, group1_indici, group2_indici, group1_label, group2_label,
                                expression = T,  differential_expression= T, use_adjusted = T, normalize_quantiles = F){
  
  #Retrieve settings
  evaluated_args <- c(as.list(environment()))
  settings <- as.list(stackoverflow::match.call.defaults()[-1])
  replace_args <- names(settings)[!names(settings) %in% unevaluated_args]
  for (argument in replace_args) {
    settings[[which(names(settings) == argument)]] <- evaluated_args[[which(names(evaluated_args) == 
                                                                              argument)]]
  }
  #Initialize outputs
  diff_genes <- NULL
  collapsed_exprs_mat <- NULL
  edgeR_deg_table <- NULL
  #Put here quantile normalization and variance stability
  if (expression){
    collapsed_exprs_mat <- DESeq2::varianceStabilizingTransformation(count_matrix)
    if(normalize_quantiles){
      collapsed_exprs_mat <- preprocessCore::normalize.quantiles(as.matrix(collapsed_exprs_mat), copy = TRUE)
    }
  }
  #Same for expression data
  if (differential_expression == T){
    P.Value = NULL
    adj.P.Val = NULL
    
    group_factor <- create_group_factor(samples = colnames(count_matrix),
                                        group1_indici = group1_indici,
                                        group2_indici = group2_indici)
    
    y <- edgeR::DGEList(counts = count_matrix, group = group_factor)
    
    y <- edgeR::calcNormFactors(y)
    
    design <- model.matrix(~group_factor)
    
    y <- edgeR::estimateDisp(y, design, robust = T)
    
    fit <- edgeR::glmQLFit(y, design)
    qlf <- edgeR::glmQLFTest(fit )
    edgeR_deg_table <- edgeR::topTags(qlf, n = nrow(qlf$counts))
    if (use_adjusted){
      diff_genes <- data.frame(rownames(edgeR_deg_table$table), edgeR_deg_table$table$FDR)
    }else{
      diff_genes <- data.frame(rownames(edgeR_deg_table$table), edgeR_deg_table$table$PValue)
    }
    diff_genes <- stats::na.omit(diff_genes)
    
    
  }
  modifier_input <- create_custom_rna_input_object(diff_genes = diff_genes, 
                                                   edgeR_deg_table = edgeR_deg_table$table,
                                                   annotated_exprs_matrix = collapsed_exprs_mat,
                                                   count_matrix = count_matrix,
                                                   group1_indici = group1_indici,
                                                   group2_indici = group2_indici,
                                                   group1_label = group1_label,
                                                   group2_label = group2_label,
                                                   settings = settings)
  
  return (modifier_input)
}

#Calculate differentially expressed genes for microarray data
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
#Helper function for differential expression for microarrays, adds annotation to probes
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
#' Create a generic microarray based input object
#' @inheritParams create_input_microarray
#' @param diff_genes A 2 two column data.frame where the first column are genes and the second column are p-values
#' @param limma_probe_table A data.frame from \code{limma topTable} with added gene annotation
#' @param annotated_exprs_matrix A matrix where the rows are genes and the columns samples.
#' @param settings Settings used to generate the object. Used only internally by the package.
#' 
#' @details 
#' This function allows the creation of a generic microarray based input object with the same class as objects created by
#' \code{\link{create_input_microarray}}. This can be useful in the cases where you already have differentially 
#' expressed genes, an annotated expression matrix or both and want to wrap that into an input object
#' to use in downstream analysis. All arguments are optional.
#' @seealso 
#' \code{\link{create_input}}
#' @return
#' The function returns an object of class "MODifieR_input". The object is a named list containing the
#' following components:
#' \item{diff_genes}{A 2 two column data.frame where the first column are genes and the second column p-values}
#' \item{limma_probe_table}{A data.frame from \code{limma topTable} with added gene annotation}
#' \item{annotated_exprs_matrix}{A matrix where the rows are genes and the columns samples. Probes have been collapsed
#' into genes using \code{collapse_method}}
#' \item{expression_matrix}{A matrix, the original input expression matrix}
#' \item{annotation_table}{A data.frame, the original annotation table used to annotate the probes}
#' \item{group_indici}{A named list containing 2 numeric vectors. The names are the group labels and the values 
#' are the group indici}
#' 
#' @author Dirk de Weerd
#' 
#' @export
create_custom_microarray_input_object <- function(diff_genes = NULL, limma_probe_table = NULL,
                                                  annotated_exprs_matrix = NULL, expression_matrix = NULL, 
                                                  annotation_table = NULL, group1_indici = NULL,
                                                  group2_indici = NULL, group1_label = NULL, 
                                                  group2_label = NULL, settings = NULL){
  evaluated_args <- c(as.list(environment()))
  validated_input <- do.call("validate_new_input_objects", evaluated_args)
  
  group_indici <- list("group_1_indici" = validated_input$group1_indici,
                       "group_2_indici" = validated_input$group2_indici)
  names(group_indici) <- c(validated_input$group1_label, validated_input$group2_label)
  
  modifier_input <- list("diff_genes" = validated_input$diff_genes,
                         "limma_probe_table" = limma_probe_table,
                         "annotated_exprs_matrix" = validated_input$annotated_exprs_matrix,
                         "expression_matrix" = expression_matrix,
                         "annotation_table" = annotation_table,
                         "group_indici" = group_indici,
                         "settings" = settings)
  class(modifier_input) <- c("MODifieR_input", "Expression", "MicroArray")
  
  return (modifier_input)
  
  
}
#
summarize_probes <- function(limma_probe_table, use_adjusted){
  adj.P.Val <- NULL
  P.Value <- NULL
  if (use_adjusted){
    diff_genes <- plyr::ddply(.data = , limma_probe_table,
                              .variables = "IDENTIFIER", .fun = plyr::summarise, pvalue = min(adj.P.Val))
  }else{
    diff_genes <- plyr::ddply(.data = , limma_probe_table,
                              .variables = "IDENTIFIER", .fun = plyr::summarise, pvalue = min(P.Value))
  }
  return (diff_genes)
}

collapse_probes <- function(expression_matrix, annotation_table, method){
  
  freq_probes <- table(annotation_table$PROBEID)
  
  probe_names <- names(freq_probes[freq_probes == 1])
  
  annotation_table_expression <- annotation_table[annotation_table$PROBEID  %in% probe_names,]
  
  annotation_table_expression <- stats::na.omit(annotation_table_expression)
  
  collapsed_data <- WGCNA::collapseRows(datET = expression_matrix,
                                        rowGroup =  annotation_table_expression$IDENTIFIER,
                                        rowID =  annotation_table_expression$PROBEID, method = method)
  
  collapsed_exprs_mat <- collapsed_data$datETcollapsed
  
  return (collapsed_exprs_mat)
}
#' Recalculate DEGs 
#' @inheritParams create_input_microarray
#' @inheritParams clique_sum_permutation
#' 
#' @details 
#' Recalculate DEGs to either use adjusted or unadjusted p values 
#' 
#' @return 
#' 
#' MODifieR_input object
#' 
#' @seealso
#' 
#' \code{\link{create_input}}
#' 
#' @author Dirk de Weerd
#' 
#' @export
recalculate_diff_genes <- function(MODifieR_input, use_adjusted){
  MODifieR_input$diff_genes <- summarize_probes(limma_probe_table = MODifieR_input$limma_probe_table, 
                                                use_adjusted = use_adjusted)
  MODifieR_input$settings$use_adjusted <- use_adjusted
  
  return(MODifieR_input)
}
#' Recalculate collapsing probes to genes
#' @inheritParams create_input_microarray
#' @inheritParams clique_sum_permutation
#' 
#' @details 
#' Recalculate the collapsing of probes to genes using on the \code{method} options
#' 
#' @return 
#' 
#' MODifieR_input object
#' 
#' @seealso
#' 
#' \code{\link{create_input}}
#' 
#' @author Dirk de Weerd
#' 
#' @export
recalculate_expression <- function(MODifieR_input, method){
  MODifieR_input$annotated_exprs_matrix <- collapse_probes(expression_matrix = MODifieR_input$expression_matrix, 
                                                           annotation_table = MODifieR_input$annotation_table, 
                                                           method = method)
  MODifieR_input$settings$method <- method
  MODifieR_input$settings$expression <- TRUE
  return(MODifieR_input)
}
#' Create a generic RNA-seq based input object
#' @inheritParams create_input_microarray
#' @inheritParams create_input_rnaseq
#' @param diff_genes A 2 two column data.frame where the first column are genes and the second column are p-values
#' @param edgeR_deg_table A data.frame from \code{edgeR glmQLFit}
#' @param annotated_exprs_matrix A matrix where the rows are genes and the columns samples, used for 
#' WGCNA-based methods
#' @param settings Settings used to generate the object. Used only internally by the package.
#' 
#' @details 
#' This function allows the creation of a generic RNA-seq based input object with the same class as objects created by
#' \code{\link{create_input_rna}}. This can be useful in the cases where you already have differentially 
#' expressed genes, an annotated expression matrix or both and want to wrap that into an input object
#' to use in downstream analysis. All arguments are optional.
#' @seealso 
#' \code{\link{create_input_rnaseq}}
#' @return
#' The function returns an object of class "MODifieR_input". The object is a named list containing the
#' following components, given that they are provided as arguments to function first:
#' \item{diff_genes}{A 2 two column data.frame where the first column are genes and the second column unadjusted p-values}
#' \item{edgeR_deg_table}{A data.frame from \code{edgeR glmQLFit}}
#' \item{annotated_exprs_matrix}{A matrix where the rows are genes and the columns samples. Probes have been collapsed
#' into genes using \code{collapse_method}}
#' \item{count_matrix}{A matrix, the original input expression matrix}
#' \item{group_indici}{A named list containing 2 numeric vectors. The names are the group labels and the values 
#' are the group indici}
#' 
#' @author Dirk de Weerd
#' 
#' @export
create_custom_rna_input_object <- function(diff_genes = NULL, edgeR_deg_table = NULL,
                                           annotated_exprs_matrix = NULL, count_matrix = NULL, 
                                           group1_indici = NULL, group2_indici = NULL, group1_label = NULL, 
                                           group2_label = NULL, settings = NULL){
  
  evaluated_args <- c(as.list(environment()))
  validated_input <- do.call("validate_new_input_objects", evaluated_args)
  
  group_indici <- list("group_1_indici" = validated_input$group1_indici,
                       "group_2_indici" = validated_input$group2_indici)
  names(group_indici) <- c(validated_input$group1_label, validated_input$group2_label)
  
  modifier_input <- list("diff_genes" = validated_input$diff_genes,
                         "edgeR_deg_table" = edgeR_deg_table,
                         "annotated_exprs_matrix" = validated_input$annotated_exprs_matrix,
                         "count_matrix" = count_matrix,
                         "group_indici" = group_indici,
                         "settings" = settings)
  class(modifier_input) <- c("MODifieR_input", "Expression", "RNA-seq")
  
  return (modifier_input)
  
  
}