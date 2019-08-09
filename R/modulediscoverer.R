#'Module Discoverer
#'
#' A clique based algorithm by Vlaic et al. to produce disease module from Differentially Expressed Genes
#'
#'
#' @inheritParams clique_sum_exact
#' @inheritParams build_clique_db
#' @param repeats Number of times the algorithm is repeated
#' @param permutations  Number of permutations to perform to identify the community structure
#' @param clique_cutoff cutoff pvalue for significant cliques
#' 
#' @details 
#' This is an implementation of the \emph{single seed} Module Discoverer algorithm.
#' The code has been adapted from the orignal code by Vlaic et al. For details, please see the paper referenced below
#' 
#' @return 
#' modulediscoverer returns an object of class "MODifieR_module" with subclass "module_discoverer". 
#' This object is a named list containing the following components:
#' \item{module_genes}{A character vector containing the genes in the final module}
#' \item{graph}{\code{\link[igraph]{igraph}} graph containing the disease module}
#' \item{settings}{A named list containing the parameters used in generating the object}
#' 
#' @seealso 
#' 
#' \url{https://www.leibniz-hki.de/en/modulediscoverer.html}
#' 
#' @references \cite{Vlaic, S., Tokarski-schnelle, C., Gustafsson, M., Dahmen, U., Guthke, R., & Schuster, S. (2017). 
#' ModuleDiscoverer: Identification of regulatory modules in protein-protein interaction networks., 1–17.}
#' @export
modulediscoverer <- function(MODifieR_input, ppi_network, permutations = 10000, deg_cutoff = 0.05, repeats = 15,
                             clique_cutoff = 0.01, dataset_name = NULL){
  
  # Retrieve settings
  evaluated_args <- c(as.list(environment()))
  settings <- as.list(stackoverflow::match.call.defaults()[-1])
  replace_args <- names(settings)[!names(settings) %in% unevaluated_args]
  for (argument in replace_args) {
    settings[[which(names(settings) == argument)]] <- evaluated_args[[which(names(evaluated_args) == 
                                                                              argument)]]
  }
  #Validate the input parameters
  check_diff_genes(MODifieR_input, deg_cutoff = deg_cutoff)
  ppi_network <- validate_ppi(ppi_network)
  validate_inputs(settings)
  
  if (!is.null(dataset_name)){
    settings$MODifieR_input <- dataset_name
  }
  
  diffgen <- MODifieR_input$diff_genes
  colnames(diffgen) <- c("gene" , "p_val")
  ppi_network <- graph.data.frame(ppi_network , directed = FALSE)
  ppi_network <- simplify(ppi_network, remove.multiple = TRUE, remove.loops = TRUE)
  
  A <- get.adjacency(ppi_network)
  
  proteins <- rownames(A)
  colnames(A) <- NULL
  rownames(A) <- NULL
  
  degrees <- degree(ppi_network)
  
  vlist <- cbind("content" = 1:nrow(A), "weight" = rep(1, nrow(A)), "degree" = degrees)
  
  background <- unique(diffgen[,1])
  
  degs <- diffgen$gene[diffgen$p_val < deg_cutoff]
  
  degs_random_datasets <- replicate(n = permutations, expr = sample(background,
                                                                    size = length(degs),
                                                                    replace = FALSE), simplify = F)
  
  cl = parallel::makeCluster(3) # initialize the cluster with x cores.
  doParallel::registerDoParallel(cl)
  
  parallel::clusterCall(cl, function(x) .libPaths(x), .libPaths())
  
  db_results_singleSeed <- foreach(j = 1:repeats, .combine = 'append', .packages = "MODifieRDev") %dopar% {
    cat(paste("processing run:",j,'\n'))
    set.seed(j) # if we don't set a seed, each repeat will return identical results.
    db_results <- lapply(1:permutations, function(i){return(MODifieRDev:::moduleDiscoverer.fragmentGraph(A=A, vlist=vlist, nbrOfSeeds=1))})
    return(db_results)
  }
  
  parallel::stopCluster(cl) # stop the cluster
  
  database_singleSeed <- moduleDiscoverer.createDatabase(results=db_results_singleSeed,
                                                         proteins=proteins)
  
  input_singleSeed = moduleDiscoverer.db.create_MD_object(database=database_singleSeed,
                                                          foregrounds=list("NASH"=degs),
                                                          cores=2, background=background,
                                                          chunks=100,
                                                          randomDataSets=list(degs_random_datasets))
  
  result_singleSeed = moduleDiscoverer.db.testForCliqueEnrichment(database=database_singleSeed,
                                                                  input=input_singleSeed, cores=3)
  
  result_singleSeed.ec = moduleDiscoverer.db.extractEnrichedCliques(database=database_singleSeed,
                                                                    result=result_singleSeed,
                                                                    p.value=clique_cutoff)
  
  module_singleSeed = moduleDiscoverer.module.createModule(result=result_singleSeed.ec)
  
  if (!is.null(module_singleSeed)){
    vertices <- names(V(module_singleSeed))
  }else{
    vertices <- NULL
  }
  
  new_modulediscoverer_module <- list("module_genes" =  vertices,
                                      "graph" = module_singleSeed,
                                      "settings" = settings)
  
  class(new_modulediscoverer_module) <- c("MODifieR_module", "module_discoverer", class(MODifieR_input)[3])
  
  return(new_modulediscoverer_module)
  
}
