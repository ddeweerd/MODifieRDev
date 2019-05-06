#' Clique Sum Exact
#'
#' An implementation of the clique-based disease module inference method proposed by Barrenäs et al. 
#' @inheritParams clique_sum_permutation
#' @param MODifieR_input A MODifieR input object produced by \code{\link{create_input}} function
#' @param db A clique database created by \code{\link{build_clique_db}} 
#' @param clique_significance p-value for cliques to be considered significant
#' @param deg_cutoff  p-value cutoff for differentialy expressed genes
#' @param min_clique_size Minimal size for cliques 
#' @param min_deg_in_clique Minimum number of DEGs to be present in a clique
#' @param multiple_cores Parallel process using multiple cores?
#' @param n_cores Number of cores to use if parallel processing
#' @param dataset_name Optional name for the input object that will be stored in the settings object.
#' Default is the variable name of the input object
#' @return clique_sum_exact returns an object of class "MODifieR_module" with subclass "Clique_Sum_exact". 
#' This object is a named list containing the following components:
#' \item{module_genes}{A character vector containing the genes in the final module}
#' \item{settings}{A named list containing the parameters used in generating the object}
#' @details 
#'  Clique_sum_exact finds cliques of at least size \code{min_clique_size} that are 
#'  significantly enriched with DEGs. 
#'  The union of maximal cliques with a Fisher-exact test p-value below \code{clique_significance} and at 
#'  least \code{min_deg_in_clique} is the final disease module.
#' @references 
#' \cite{Barrenäs F, Chavali S, Alves AC, et al. Highly interconnected genes in 
#' disease-specific networks are enriched for disease-associated polymorphisms. 
#' Genome Biology. 2012;13(6):R46. doi:10.1186/gb-2012-13-6-r46.}
#' @export
clique_sum_exact <- function(MODifieR_input, db, clique_significance = 0.01,
                             deg_cutoff = 0.05, min_clique_size = 5, min_deg_in_clique = 3, 
                             multiple_cores = T, n_cores = 4, dataset_name = NULL){
  
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
  
  deg_genes <- genes[genes > -log10(deg_cutoff)]
  non_deg_genes <- genes[genes < -log10(deg_cutoff)]
  
  #Fix tables
  tables <- RSQLite::dbListTables(con)
  
  chunk_table <- tables[grep(pattern = "chunk_summary", x = tables)]
  
  tables <- tables[-grep(pattern = "chunk_summary", x = tables)]
  tables <- tables[-grep(pattern = "unique_genes", x = tables)]
  
  table_numbers <- as.numeric(gsub(pattern = "clique", replacement = "", x = tables))
  
  tables <- tables[as.numeric(gsub(pattern = "clique", replacement = "", x = tables)) >= min_clique_size]
  
  tables <- tables[order(as.numeric(gsub(pattern = "clique", replacement = "", x = tables)))]
  
  table_sizes <- as.numeric(sub("clique", "", tables))
  
  cutoff_per_size <- get_cutoffs(min_clique_size = min_clique_size, 
                                 max_table_size = max(table_sizes), 
                                 total_n_genes = length(genes), 
                                 deg_genes = deg_genes, 
                                 unsig = length(non_deg_genes), 
                                 clique_significance = clique_significance, 
                                 min_deg_in_clique = min_deg_in_clique)
  
  chunk_table <- get_chunk_table(con = con, 
                                 deg_genes = deg_genes, 
                                 cutoff_per_size = cutoff_per_size,
                                 min_clique_size = min_clique_size,
                                 genes = genes,
                                 min_deg_in_clique = min_deg_in_clique)
  
  if (is.null(chunk_table)){
    module_genes <- NULL
  }else{
    if (multiple_cores == T){
      cl <- parallel::makeCluster(n_cores, outfile = "")
      doParallel::registerDoParallel(cl)
      parallel::clusterCall(cl, function(x) .libPaths(x), .libPaths())
      module_genes <- foreach(i = 1:length(chunk_table),.combine = 'append', .packages = "MODifieRDev" ) %dopar% {
        clique_sum_core_exact(chunk_row = chunk_table[[i]], db = db, deg_genes = deg_genes, 
                              cutoff_per_size = cutoff_per_size, 
                              genes = genes,
                              min_clique_size = min_clique_size)
      }
        module_genes <- unique(unlist(module_genes))
        parallel::stopCluster(cl) # stop the cluster
      
    }else{
      #Convert list to table
      chunk_table <- do.call(rbind, chunk_table)
      
      module_genes <- unique(unlist(apply(X = chunk_table, 1, clique_sum_core_exact,  
                                          db = db, deg_genes = deg_genes, genes = genes, cutoff_per_size = cutoff_per_size, min_clique_size = min_clique_size)))
    }
  }
  
  #Construct list and give it the correct class
  new_clique_sum_module <- list("module_genes" =  module_genes,
                                "settings" = settings)
  
  class(new_clique_sum_module) <- c("MODifieR_module", "Clique_Sum_exact")
  
  RSQLite::dbDisconnect(conn = con)
  
  return(new_clique_sum_module)
}
#Create a set of 4 column tables for all clique sizes.
#Column 1: Number of DEGs in clique
#Column 2: Number of DEGs not in clique
#Column 3: Number of unsignificant genes in clique
#Column 4: Number of unsignificant genes not in clique
#The tables are used to calculate the minimum number of DEGs needed
#for every clique size to achieve significance
create_fisher_tables <- function(clique_size, total_n_genes, n_sig, unsig){
  
  min_size <- min(c(clique_size, n_sig))
  fisher_table <- cbind(1:min_size, clique_size - 1:min_size, n_sig - 1:min_size, unsig - clique_size + 1:min_size)
  
  fisher_table[which(fisher_table < 0)] <- 0
  
  return(fisher_table)
}
#Evaluate a 'Fisher table' and return the minimum number of DEGs needed to achieve significance.
#If this is less than min_deg_in_clique, return min_deg_in_clique
evaluate_table <- function(fisher_table, clique_significance, min_deg_in_clique){
  for (i in 1:nrow(fisher_table)){
    if (MODifieRDev:::fisher_pval(clique_row = fisher_table[i,]) <= clique_significance)
      return(max(i, min_deg_in_clique))
  }
    return (100000)
  
}
#Function that determines if there are potentially significant genes in a clique
get_chunks <- function(result_row, cutoff_per_size, deg_genes, genes, min_clique_size, min_deg_in_clique){
  #Initialize variable 'cutoff' for current chunk
  deg_cutoff <- NULL
  #Get the row genes
  row_genes <- unname(unlist(strsplit(x = as.character(result_row[4]), split = " ")))
  #Get the number of genes in the chunk that are actually present in the data
  n_p_genes <-  sum(names(genes) %in% row_genes)
  #Quantifying the number of missing genes, i.e. genes that are present in the chunk, but 
  #not in the data
  n_missing_genes <- length(row_genes) - n_p_genes 
  #number of missing genes will say something about the potential clique size. For example, 
  #if clique size is 65, and there are 55 genes missing, that means there are potential size 10 cliques.
  smallest_clique <- unlist(result_row[3] - n_missing_genes)
  #If the number of missing genes is equal or more than the clique size, assumption is that there are
  #potential cliques of size min_clique_size cliques DEG cutoffs should be set accordingly
  if (smallest_clique <= min_clique_size){
    deg_cutoff <- cutoff_per_size[[1]]
    #If the number of missing genes is subtracted from the clique size and there is still a minimum clique left
    #that is bigger than min_clique_size, assumption is that there is a potential clique of size
    # clique_size - missing genes. DEG cutoffs should be set accordingly
  }else{
    deg_cutoff <- cutoff_per_size[as.character(smallest_clique)]
  }
  #Count the number of DEGs present in the chunk
  n_deg_genes <- sum(names(deg_genes) %in% row_genes)
  #If the number of DEGs in the chunk is larger than the cutoff, there are potentially significant cliques.
  #Return the clique_size (that is the table in the DB), the start and stop of the chunk. 
  #Later this can be converted in a query. If the number of DEGs is lower than the cutoff, do nothing.
  if (n_deg_genes >= deg_cutoff){
    return (c(paste("clique", result_row[3], sep = ""), result_row[1], result_row[2], deg_cutoff))
  }
}

get_chunk_table <- function(con, deg_genes, cutoff_per_size, min_clique_size, genes, min_deg_in_clique){
  #Select all chunk summaries where minimum clique size is at least min_clique_size
  query <- sprintf("SELECT * FROM chunk_summary WHERE clique_size >= %s ORDER BY clique_size", min_clique_size)
  result_table <- RSQLite::dbGetQuery(con, query)
  #Get a list of chunks with potentially significant genes
  chunk_list <- plyr::compact(plyr::alply(result_table, 1, get_chunks, 
                                          cutoff_per_size = cutoff_per_size, 
                                          deg_genes = deg_genes, 
                                          genes = genes, 
                                          min_clique_size = min_clique_size, 
                                          min_deg_in_clique = min_deg_in_clique))
  
  
  return (chunk_list)
  
}
#Core algorithm, recieves a 'chunk_row' vector containing a:
#clique table [1]
#start coordinate [2]
#stop coordinate [3]
#A database query is constructed with these values. Executing the query will give a set of cliques
#that will be checked for significance. Return a unique set of genes from the siginificant cliques
clique_sum_core_exact <- function(chunk_row, db, deg_genes, genes, cutoff_per_size, min_clique_size){
  
  con <- RSQLite::dbConnect(RSQLite::SQLite(), db)
  query <- sprintf("SELECT * FROM %s WHERE rowid BETWEEN %s AND %s", chunk_row[1], chunk_row[2], chunk_row[3])
  result <- RSQLite::dbGetQuery(con, query) 
  #Make a 2 column matrix with:
  #column 1: Number of genes present in the input dataset
  #column 2: Number of DEGs in these genes
  length_n_degs <- cbind(apply(X = result, MARGIN = 1, FUN = function(x)sum(as.character(x) %in% names(genes))), 
                         apply(X = result, MARGIN = 1, FUN = function(x)sum(as.character(x) %in% names(deg_genes))))
  #Check for each row if the clique is significant. Return type is a logical vector
  significant_rows <- apply(X = length_n_degs, MARGIN = 1, FUN = get_significant_rows, 
                            cutoff_per_size = cutoff_per_size, min_clique_size = min_clique_size)
  #Subset the result with the logical vector
  module_genes <- result[significant_rows, ]
  #Return the unique genes in vector format
  
  RSQLite::dbDisconnect(conn = con)
  
  
  return (unique(unlist(module_genes)))
  
}

get_significant_rows <- function(length_deg_row, cutoff_per_size, min_clique_size){
  #If the number of genes in the clique is the same or higher than min_clique_size
  if (length_deg_row[1] >= min_clique_size){
    #If the number of DEGs is sufficient to make clique significant return TRUE
    if (length_deg_row[2] >= cutoff_per_size[as.character(length_deg_row[1])]){
      return(T)
      #Else return F
    }else{
      return(F)
    }
  }else{
    return(F)
  }
}
#Retrieve for every possible clique size (min_clique_size to size of largest clique) the number of genes that is needed
#to get a p-value under clique_significance
get_cutoffs <- function(min_clique_size, max_table_size, total_n_genes, deg_genes, unsig, clique_significance, min_deg_in_clique){
  fisher_tables <- lapply(X = min_clique_size:max_table_size, create_fisher_tables, 
                          total_n_genes = total_n_genes, 
                          n_sig = length(deg_genes), 
                          unsig = unsig)
  #The names of the list are the sizes of the cliques; the values are the number of DEGs needed to acieve significance
  cutoff_per_size <- lapply(X = fisher_tables, FUN = evaluate_table, clique_significance = clique_significance, min_deg_in_clique)
  names(cutoff_per_size) <- paste0(min_clique_size:max_table_size)
  #If the size is NULL there is no number of DEGs possible to achieve significance, set to arbitrary high number instead
  cutoff_per_size[sapply(cutoff_per_size, is.null)] <- 10000
 
  return(unlist(cutoff_per_size))
}