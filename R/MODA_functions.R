#Adapted code from original MODA package
WMPH <- function (datExpr, indicatename, cutmethod = c("Density", "Modularity"), power = 10) 
{
  #dir.create(file.path("./", foldername), showWarnings = FALSE)
  ADJ1 = abs(cor(datExpr, use = "p"))^power
  dissADJ = 1 - ADJ1
  hierADJ = flashClust::hclust(stats::as.dist(dissADJ), method = "average")
  cutHeights <- seq(0.05, 1, by = 0.05)
  NumCutHeights <- length(cutHeights)
  pDensity <- numeric(length <- NumCutHeights)
  maxpDensity <- 0
  maxHeight <- 0
  for (i in 1:NumCutHeights) {
    groups <- stats::cutree(hierADJ, h = cutHeights[i])
    if (cutmethod == "Density") {
      pDensity[i] <- MODA::PartitionDensity(ADJ1, groups)
    }
    else {
      pDensity[i] <- MODA::PartitionModularity(ADJ1, groups)
    }
    if (pDensity[i] > maxpDensity) {
      maxpDensity <- pDensity[i]
      maxHeight <- cutHeights[i]
    }
  }
  minModuleSize = 30
  dynamicMods = dynamicTreeCut::cutreeDynamic(dendro = hierADJ, distM = dissADJ, 
                              cutHeight = maxHeight, deepSplit = 2, pamRespectsDendro = FALSE, 
                              minClusterSize = minModuleSize)
  dynamicColors = WGCNA::labels2colors(dynamicMods)
  intModules = table(dynamicColors)

  densegenes_list <- list()
  for (J in 1:length(intModules)) {
    #idx is a vector of colnumbers (matrix is transposed), 
    #corresponding to the genes in DenseGenes
    #intModules is iterated over, into idx which is in turn used to subset, so
    #number of modules = length of densegenes_list
    idx <- which(dynamicColors == names(intModules)[J])
    DenseGenes = colnames(datExpr)[idx]
  
    densegenes_list[[J]] <- DenseGenes
  }
  
  
  
  return(densegenes_list)
}
#Adapted code from original MODA package
compare_modules <- function (module_list1, module_list2) 
{
  #Set the lengths for simplicity
  nm1 <- length(module_list1)
  nm2 <- length(module_list2)
  #Two dimensional array according to size of modules, first value (nm1) is number of rows
  #Second value, nm2, is number of cols.
  my.array <- array(0, dim = c(nm1, nm2))
  for (i1 in 1:nm1) {
    for (i2 in 1:nm2) {
        my.array[i1, i2] = length(intersect(module_list1[[i1]], module_list2[[i2]])) / 
          length(union(module_list1[[i1]], module_list2[[i2]]))
    }
  }
  return(my.array)
}
#Adapted code from original MODA package
moda_extract_modules_index_specific <- function(jaccard_table, specificTheta = 0.1){
  which(rowSums(jaccard_table) <= min(rowSums(jaccard_table)) + specificTheta)
}
#Adapted code from original MODA package
moda_extract_modules_index_conserved <- function(jaccard_table, conservedTheta = 0.1){
  which(rowSums(jaccard_table) >= max(rowSums(jaccard_table)) - conservedTheta)
}
