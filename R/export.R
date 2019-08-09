#' Export MODifieR objects to csv files
#' 
#' Writes all the fields of a MODifieR object to separate csv files
#' 
#' @param MODifieR_object The object to export
#' @param folder The folder to write to
#' 
#' @export
export <- function(MODifieR_object, folder){
  MODifieR_object <- prepare_object(MODifieR_object = MODifieR_object)
  if (length(grep(pattern = "/$", folder)) == 0){
    folder <- paste0(folder, "/")
  }
  MODifieR_object <- prepare_export(MODifieR_object = MODifieR_object)
  write_export(MODifieR_object = MODifieR_object, folder = folder)
}
#' Export MODifieR to object to xlsx format
#' 
#' Writes a MODifieR object to a xlsx workbook with multiple sheets
#' @inheritParams export
#' @param filename The name of xlsx file tow rite to. If NULL, the name of the object will be used
#' @export
export_xlsx <- function(MODifieR_object, folder, filename = NULL){
  
  dir.create(folder, showWarnings = F, recursive = T)
  
  MODifieR_object <- prepare_object(MODifieR_object = MODifieR_object)
  
  args <- as.list(stackoverflow::match.call.defaults()[-1])
  if (length(grep(pattern = "/$", folder)) == 0){
    folder <- paste0(folder, "/")
  }
  
  if (is.null(filename)){
    filename  <- args$MODifieR_object
  }
  
  xlsx_file <- paste0(folder, filename, ".xlsx")
  
  MODifieR_object <- preprocess_object(MODifieR_object)
  class(MODifieR_object) <- "list"
  write.xlsx(x = MODifieR_object, file = xlsx_file)
}

write_export <- function(MODifieR_object, folder){
  dir.create(folder, showWarnings = F, recursive = T)
  for (i in which(names(MODifieR_object) %in% remove_headers)){
    write.table(x = MODifieR_object[[i]], file = paste0(folder, names(MODifieR_object)[i], ".csv"), 
                row.names = F, col.names = F, quote = F)
  }
  for (i in which(!names(MODifieR_object) %in% remove_headers)){
    write.table(x = MODifieR_object[[i]], file = paste0(folder, names(MODifieR_object)[i], ".csv"), 
                row.names = F, quote = F)
  }
}



