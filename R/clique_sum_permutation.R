#' Clique Sum
#'
#' An implementation of the clique-based disease module inference method proposed by Gustafsson et al. 
#'
#' @inheritParams clique_sum_exact
#' @param MODifieR_input A MODifieR input object produced by \code{\link{create_input}} function
#' @param n_iterations Number of iterations to be performed for the permutation based p-value
#' @param clique_significance p-value for cliques to be considered significant
#' @param min_clique_size Minimal size for cliques 
#' @param dataset_name Optional name for the input object that will be stored in the settings object.
#' Default is the variable name of the input object
#' @return clique_sum_permutation returns an object of class "MODifieR_module" 
#' with subclass "Clique_Sum_permutation". 
#' This object is a named list containing the following components:
#' \item{module_genes}{A character vector containing the genes in the final module}
#' \item{settings}{A named list containing the parameters used in generating the object}
#' @details 
#'  Clique_sum_permutation finds cliques of at least size \code{min_clique_size} that are 
#'  significantly enriched with DEGs. For every clique size, a null distribution is created using
#'  the summed -log 10 p-values. The union of maximal cliques with a summed -log 10 p-value below \code{clique_significance} and at 
#'  least \code{min_deg_in_clique} is the final disease module.
#' @references 
#' \cite{Gustafsson, M., Edström, M., Gawel, D., Nestor, C. E., Wang, H., Zhang, H., … Benson, M. (2014). 
#' Integrated genomic and prospective clinical studies show the importance of modular pleiotropy for disease susceptibility, 
#' diagnosis and treatment. Genome Medicine, 6(2), 17. https://doi.org/10.1186/gm534}
#' @export
clique_sum_permutation <- function(MODifieR_input, db, n_iterations = 10000, clique_significance = 0.01, 
                                   min_clique_size = 5,  multiple_cores = T, n_cores = 4, dataset_name = NULL){
  # Retrieve settings
  evaluated_args <- c(as.list(environment()))
  settings <- as.list(stackoverflow::match.call.defaults()[-1])
  replace_args <- names(settings)[!names(settings) %in% unevaluated_args]
  for (argument in replace_args) {
    settings[[which(names(settings) == argument)]] <- evaluated_args[[which(names(evaluated_args) == 
                                                                              argument)]]
  }
 
  if (!is.null(dataset_name)){
    settings$MODifieR_input <- dataset_name
  }
  con <- RSQLite::dbConnect(RSQLite::SQLite(), db)
  
  unique_genes <- unname(unlist(RSQLite::dbGetQuery(con, "SELECT * FROM unique_genes")))
  
  genes <- MODifieR_input$diff_genes
  
  
  #Convert 2 column dataframe to named vector, names will be the gene names, values will be the p value
  genes <- -log10(dataframe_to_vector(as.data.frame(genes)))
  
  genes <- genes[names(genes) %in% unique_genes]
  
  #Fix tables
  tables <- RSQLite::dbListTables(con)
  
  chunk_table <- tables[grep(pattern = "chunk_summary", x = tables)]
  
  tables <- tables[-grep(pattern = "chunk_summary", x = tables)]
  tables <- tables[-grep(pattern = "unique_genes", x = tables)]
  
  table_numbers <- as.numeric(gsub(pattern = "clique", replacement = "", x = tables))
  
  tables <- tables[as.numeric(gsub(pattern = "clique", replacement = "", x = tables)) >= min_clique_size]
  
  tables <- tables[order(as.numeric(gsub(pattern = "clique", replacement = "", x = tables)))]
  
  table_sizes <- as.numeric(sub("clique", "", tables))
  
  null_scores <- sapply(min_clique_size:max(table_sizes), permute_scores, n_iterations = n_iterations, genes = genes, clique_significance =clique_significance)
  
  names(null_scores) <- min_clique_size:max(table_sizes)
  
  chunk_table <- get_chunk_table_perm(con = con, 
                                      null_scores = null_scores, 
                                      min_clique_size = min_clique_size,
                                      genes = genes)
  
  
  if (is.null(chunk_table)){
    module_genes = NULL
  }else{
    if (multiple_cores == T){
      cl <- parallel::makeCluster(n_cores, outfile = "")
      doParallel::registerDoParallel(cl)
      parallel::clusterCall(cl, function(x) .libPaths(x), .libPaths())
      
      module_genes <- foreach(i = 1:length(chunk_table),.combine = 'append', .packages = "MODifieRDev" ) %dopar% {
        clique_sum_core_permutation(chunk_row = chunk_table[[i]], db = db, 
                                                  null_scores = null_scores, 
                                                  genes = genes, 
                                                  min_clique_size = min_clique_size)
        
      }
      parallel::stopCluster(cl) # stop the cluster
      module_genes <- unlist(unique(module_genes))
      
    }else{
      #Convert list to table
      chunk_table <- do.call(rbind, chunk_table)
      module_genes <- unique(unlist(apply(X = chunk_table, 1, clique_sum_core_permutation,  
                                          db = db, 
                                          null_scores = null_scores, 
                                          genes = genes, 
                                          min_clique_size = min_clique_size)))
    }
  }
  #Construct list and give it the correct class
  new_clique_sum_module <- list("module_genes" =  module_genes,
                                "settings" = settings)
  
  class(new_clique_sum_module) <- c("MODifieR_module", "Clique_Sum_permutation")
  
  RSQLite::dbDisconnect(conn = con)
  
  return(new_clique_sum_module)
}

clique_sum_core_permutation <- function(chunk_row, db, null_scores, genes, min_clique_size){
  con <- RSQLite::dbConnect(RSQLite::SQLite(), db)
  query <- sprintf("SELECT * FROM %s WHERE rowid BETWEEN %s AND %s", chunk_row[1], chunk_row[2], chunk_row[3])
  result <- RSQLite::dbGetQuery(con, query) 
  result <- data.frame(lapply(result, as.character), stringsAsFactors = F)
  significant_rows <- apply(X = result, MARGIN = 1, FUN = compare_clique_scores, 
                            genes = genes, null_scores = null_scores, min_clique_size = min_clique_size)
  
  
  module_genes <- result[significant_rows, ]
  #Return the unique genes in vector format
  RSQLite::dbDisconnect(con)
  return (unique(unlist(module_genes)))
  
}
permute_scores <- function(clique_size, n_iterations, genes, clique_significance){
  
  permutations <- replicate(n = n_iterations, expr = sum(genes[sample(length(genes), clique_size, replace = F)]))
  permutations <- sort(permutations, decreasing = T)
  cutoff_index <- round(x = n_iterations * clique_significance, digits = 0)
  cutoff_clique_size <- permutations[cutoff_index]
  return(cutoff_clique_size)
}

#Function that determines if there are potentially significant genes in a clique
get_chunks_perm <- function(result_row, null_scores, genes, min_clique_size){
  #Initialize variable 'cutoff' for current chunk
  deg_cutoff <- NULL
  
  sorted_genes <- sort(genes, decreasing = T)
  #Get the row genes
  row_genes <- unname(unlist(strsplit(x = as.character(result_row[4]), split = " ")))
  #Get the number of genes in the chunk that are actually present in the data
  n_p_genes <-  genes[names(genes) %in% row_genes]
  n_p_genes <- sort(n_p_genes, decreasing = T)
  #n_p_genes <-  intersect(x = names(genes), y =  row_genes)
  
  #Quantifying the number of missing genes, i.e. genes that are present in the chunk, but 
  #not in the data
  n_missing_genes <- length(row_genes) - length(n_p_genes) 
  #number of missing genes will say something about the potential clique size. For example, 
  #if clique size is 65, and there are 55 genes missing, that means there are potential size 10 cliques.
  smallest_clique <- unlist(result_row[3] - n_missing_genes)
  
  #If the number of missing genes is equal or more than the clique size, assumption is that there are
  #potential cliques of size min_clique_size cliques DEG cutoffs should be set accordingly
  clique_sizes <- max(c(smallest_clique, min_clique_size)):min(c(length(n_p_genes), as.numeric(result_row[3])))
  
  for (clique in clique_sizes){
    if (sum(n_p_genes[1:clique]) >= null_scores[as.character(clique)]){
      return (c(paste("clique", result_row[3], sep = ""), result_row[1], result_row[2]))
    }
  }
}

get_chunk_table_perm <- function(con, null_scores, min_clique_size, genes){
  #Select all chunk summaries where minimum clique size is at least min_clique_size
  query <- sprintf("SELECT * FROM chunk_summary WHERE clique_size >= %s ORDER BY clique_size", min_clique_size)
  result_table <- RSQLite::dbGetQuery(con, query)
  #Get a list of chunks with potentially significant genes
  chunk_list <- plyr::compact(plyr::alply(result_table, 1, get_chunks_perm, 
                                          null_scores = null_scores, 
                                          genes = genes, 
                                          min_clique_size = min_clique_size))
  
  
  return (chunk_list)
}

compare_clique_scores <- function(clique, genes, null_scores, min_clique_size){
  
  clique_genes <- genes[as.character(clique)]
  if (length(clique_genes) <= min_clique_size){
    return(F)
  }else{
    if (sum(clique_genes, na.rm = T) >= null_scores[as.character(length(clique_genes))]){
      return(T)
    }else{
      return(F)
    }
  }
}