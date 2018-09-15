process_module <- function(module){
  UseMethod("process_module", module)
}
process_module.list <- function(module){
  unlisted_module <- module[sapply(X = module, FUN = unlist_module)]
}
process_module.Clique_Sum <- function(module){
  module
}
process_module.Correlation_clique <- function(module){
  module
}
process_module.DIAMOnD <- function(module){
  module
}
process_module.module_discoverer <- function(module){
  module
}
process_module.DiME <- function(module){
  module
}
process_module.MODA <- function(module){
  module
}
process_module.WGCNA <- function(module){
  module$module_genes <- rownames(module$probe_info_table)
  return (module)
}
unlist_module <- function(module){
  UseMethod("unlist_module", module)
}
unlist_module.Mcode <- function(module){
  length(module$module_genes) > 20
}
unlist_module.DiffCoEx <- function(module){
  length(module$module_genes) > 20
}
#' @export
get_module_list <- function(result_list){
  
  maak <- lapply(result_list, process_module)
  
  module_list <- list()
  for (module in maak){
    if (!class(module)[1] == "MODifieR_module"){
      for (i in 1:length(module)){
        module_list[[length(module_list)+1]] <- module[[i]]
        
      }
    }else{
      module_list[[length(module_list)+1]] <- module
    }
  }
  
  names(module_list) <- sapply(X = module_list, FUN = MODifieRDev::extract_module_class)
  
  tabled_names <- table(names(module_list))
  
  duplicated_names <- names(tabled_names[tabled_names > 1])
  
  if (sum(tabled_names[tabled_names > 1]) != 0){
    for (i in 1:length(duplicated_names)){
      message(duplicated_names[i], " is non-unique, appending number to ", duplicated_names[i])
      names(module_list)[grep(pattern = duplicated_names[i], x = names(module_list))] <- paste0(names(module_list)[grep(pattern = duplicated_names[i], x = names(module_list))], 1:tabled_names[tabled_names > 1][i])
    }
  }
  return(module_list)
}
