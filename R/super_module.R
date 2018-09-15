#' Robust module concensus
#'
#'A final conglomeration of all modules identified by different module methods to produce all kinds of intersect modules
#'
#' @param module_list A list object with f´genes of each module
#' @param include_single_methods A boolean value , whether to include the original modules in the final result  
#' @export
super_module <- function(module_list, include_single_methods = T){
  all_modules <-  sapply(X = module_list, FUN = extract_module_data,
                         data_field = "module_genes" )
  
  if (length(module_list) != length(na.omit(names(module_list))) || (sum(names(module_list) == "")) > 0){
    message("Modules are unnamed, using inference methods as names instead")
    
    names(all_modules) <- sapply(X = module_list, FUN = extract_module_class)
  }
  
  tabled_names <- table(names(all_modules))
  
  duplicated_names <- names(tabled_names[tabled_names > 1])
  
  if (sum(tabled_names[tabled_names > 1]) != 0){
    for (i in 1:length(duplicated_names)){
      message(duplicated_names[i], " is non-unique, appending number to ", duplicated_names[i])
      names(all_modules)[grep(pattern = duplicated_names[i], x = names(all_modules))] <- paste0(names(all_modules)[grep(pattern = duplicated_names[i], x = names(all_modules))], 1:tabled_names[tabled_names > 1][i])
    }
  }
  n_modules <- length(all_modules)
  
  module_combinations <- sapply(X = 2:n_modules,
                                FUN =  get_combos, modules = names(all_modules))
  
  modulos <- invisible(sapply(X = module_combinations, FUN = helper, all_modules = all_modules))
  
  super_module <- unlist(modulos, recursive = F)
  super_module <- unlist(super_module, recursive = F)
  
  if (include_single_methods == T){
    super_module <- c(all_modules, super_module)
  }
  class(super_module) <- c("MODifieR_super_module", "intersection")
  
  return(super_module)
}

get_combos <- function(n, modules){
  combn(x = modules, m = n)
}

get_intersect <- function(modules, all_modules){
  mods <- all_modules[modules]
  tabled_frequencies <- table(unlist(mods)) / length(mods)
  n_modules <- length(mods)
  basic_name <- paste(modules, collapse = "+")
  
  #Set iterations for intersections
  iterations <- seq(from=(1/n_modules), to = 1, by = (1/n_modules))
  
  result <- list()
  
  for (i in iterations){
    current_rowname <- paste(basic_name, round(i,2), sep = "-")
    if (i == min(iterations)){
      current_rowname <- paste(basic_name, "union", sep = "-")
    }
    if (i == max(iterations)){
        current_rowname <- paste(basic_name, "intersection", sep = "-")
    }
    result[[current_rowname]] <- names(tabled_frequencies[tabled_frequencies >= i])
  }
  return(result)
}

helper <- function(x, all_modules){
  
  freq_module <- invisible(apply(X = as.matrix(x), MARGIN = 2, FUN = get_intersect, all_modules = all_modules))
  
  return(freq_module)
}

#' @export
get_intersections <- function(module_list){
  intersection_table <- table(unlist(sapply(X = module_list, FUN = function(x)x$module_genes)))
  intersections <- list()
  for (i in 1:length(module_list)){
    intersections[[i]] <- names(intersection_table[intersection_table > i])
  }
  names(intersections) <- paste("size", 1:length(module_list), "module", sep = "_")
  return(intersections)
}
