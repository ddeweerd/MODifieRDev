#'@inheritParams clique_sum
#' @export
#' 
dime <- function(MODifieR_input, dataset_name = NULL){
  
  default_args <- formals()
  user_args <- as.list(match.call(expand.dots = T)[-1])
  settings <- c(user_args, default_args[!names(default_args) %in% names(user_args)])
  
  if (!is.null(dataset_name)){
    settings$MODifieR_input <- dataset_name
  }
  
  exprs <- MODifieR_input$annotated_exprs_matrix
  nfeat <- dim(exprs)[1]
  nsamples <- dim(exprs)[2]
  exprs <- t(exprs)
  
  cor.data <- cor.jack(exprs)
  
  upper.r <- as.vector(cor.data[upper.tri(cor.data)])
  upper.ind <- which(upper.tri(cor.data))
  sorted.r <- sort((upper.r), index.return = TRUE)
  nedges <- floor(length(upper.r)*0.001)
  nedges.upper <- floor(0.5*nedges)
  nedges.lower <- nedges - nedges.upper
  selected.ind <- c(upper.ind[sorted.r$ix[1:nedges.lower]],
                    upper.ind[sorted.r$ix[(length(upper.r)-nedges.upper+1):length(upper.r)]])
  
  # Now connect these correlation pairs to form an adjacency matrix for the network,
  # and delete isolated nodes (genes with no selected coexpression partners)
  adj <- matrix(0, nfeat, nfeat)
  rownames(adj) <- as.character(rownames(cor.data))
  colnames(adj) <- as.character(colnames(cor.data))
  adj[selected.ind] <- 1
  ind <- lower.tri(adj)
  adj[ind] <- t(adj)[ind]
  toDel <- unname(which(rowSums(adj) == 0))
  adj <- adj[-toDel,-toDel]
  
  g <- igraph::graph.adjacency(adj,mode = "undirected")
  
  adj <- get.adjacency(g)
  dimVg = length(V(g))
  popsize = 50
  
  com <- commextr(d = dimVg, 
                                ps = as.integer(popsize), 
                                matrix = as.vector(adj), 
                                res = as.integer(as.vector(rep(0,dimVg))), 
                                vb = as.integer(0))
  
  probes <- (V(g)$name[com$res!=0])
  
  annotation_table <- MODifieR_input$annotation_table
  
  module_genes <- probes
  
  module_genes <- module_genes[!is.na(module_genes)]
  
  new_dime_module <- list("module_genes" =  module_genes,
                          "settings" = settings)
  
  class(new_dime_module) <- c("MODifieR_module", "DiME")
  
  
  return(new_dime_module)
}
