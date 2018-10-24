connect_geo <- function(){
  if(!file.exists('GEOmetadb.sqlite')){
    message("No local GEO database file, fetching file")
    GEOmetadb::getSQLiteFile()
  }
  else{
    message("Timestamp GEO database: ", file.info('GEOmetadb.sqlite')$mtime)
  }

  con <- RSQLite::dbConnect(RSQLite::SQLite(),'GEOmetadb.sqlite')

  return(con)

}

retrieve_geo_metadata <- function(geo_accession, full_table = F){
  con <- connect_geo()
  if (full_table == T){
    query <- paste0('SELECT * FROM gsm WHERE gsm IN (SELECT gsm FROM gse_gsm WHERE gse = "',  geo_accession, '")')
    result_table <- RSQLite::dbGetQuery(con, query)
  }
  else{
    query <- paste0('SELECT title, gsm, source_name_ch1, characteristics_ch1 FROM gsm WHERE gsm IN (SELECT gsm FROM gse_gsm WHERE gse = "',  geo_accession, '")')
    result_table <- RSQLite::dbGetQuery(con, query)
  }
  RSQLite::dbDisconnect(con)
  return(result_table)
}

retrieve_geo_gpl <- function(geo_accession){
  con <- connect_geo()
  query <- paste0('SELECT gpl FROM gsm WHERE gsm IN (SELECT gsm FROM gse_gsm WHERE gse = "',  geo_accession, '")')
  gpl <- RSQLite::dbGetQuery(con, query)
  RSQLite::dbDisconnect(con)
  return(as.vector(as.character(unname(unique(gpl)))))
}

retrieve_gpl_info <- function(gpl_accession){
  con <- connect_geo()
  query <- paste0('SELECT * FROM gpl WHERE gpl IS"',  gpl_accession, '"')
  gpl <- RSQLite::dbGetQuery(con, query)
  RSQLite::dbDisconnect(con)
  return(gpl)
}
retrieve_gsm <- function(geo_gsm_table, columns, regular_expressions){

  gsm_numbers <- list()
  for (i in 1:length(columns)){
    gsm_numbers[[i]] <- geo_gsm_table[grep(pattern = regular_expressions[i], x = geo_gsm_table[ ,columns[i] ]), "gsm"]
  }
  gsm_numbers <- Reduce(intersect, gsm_numbers)
  return(gsm_numbers)
}

#' @export
get_string_DB_ppi <- function(version = "10", score_threshold = 900, simplify_graph = T){
  string_db <- STRINGdb$new( version=as.character(version), species=9606,
                             score_threshold=score_threshold, input_directory="" )
  ppi_network <- string_db$get_graph()
  
  if (simplify_graph == T){
    ppi_network <- igraph::simplify(ppi_network)
  }
  
  return (ppi_network)
}

#' @export
symbol_to_prot <- function(dataframe, remove_NA = T){
  prot_ids_list <-  mapIds(org.Hs.eg.db,
                           keys=dataframe[,1],
                           column="ENSEMBLPROT",
                           keytype="SYMBOL",
                           multiVals="list")
  prot_ids_unlisted <- lapply(X = 1:length(prot_ids_list), FUN = function(x){
    cbind(names(prot_ids_list[x]), unlist(prot_ids_list[x], use.names = F))
  })
  
  prot_ids_final <- base::Reduce(f = rbind, x = prot_ids_unlisted)
  if (remove_NA == T){
    prot_ids_final <- na.omit(prot_ids_final)
  }
  
  colnames(prot_ids_final) <- c("hgnc_symbol", "ensembl_peptide_id")
  
  return(prot_ids_final)
}
#' @export
prot_to_entrez <- function(dataframe){
  entrez_ids <-  na.omit(t(as.data.frame(mapIds(org.Hs.eg.db,
                                                keys= dataframe[,1],
                                                column="ENTREZID",
                                                keytype="ENSEMBLPROT",
                                                multiVals="list"))))
  return(entrez_ids[,1])
}
