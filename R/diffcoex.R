#' DiffCoEx
#' 
#' A Coexpression based Algorithm to identify disease module 
#' @inheritParams clique_sum
#' @param dataset1 A dataframe of the microarray matrix of the Control samples
#' @param dataset2 A dataframe of the microarray matrix of the Patient samples
#' @param cor_method correlation method to be used (default is "spearman")
#' @export

diffcoex <- function(MODifieR_input, value = 20, beta = 6, cor_method = "spearman",
                     cluster_method = "average", cuttree_method = "hybrid",
                     cut_height = 0.996, deep_split = 0, pam_respects = F,
                     min_cluster_size = 20, cluster_merge_cut_height = 0.2, dataset_name = NULL){

  # Retrieve settings
  default_args <- formals()
  user_args <- as.list(match.call(expand.dots = T)[-1])
  settings <- c(user_args, default_args[!names(default_args) %in% names(user_args)])
  
  if (!is.null(dataset_name)){
    settings$MODifieR_input <- dataset_name
  }
  
  dataset1 <- t(MODifieR_input$annotated_exprs_matrix[,MODifieR_input$group_indici[[1]]])
  dataset2 <- t(MODifieR_input$annotated_exprs_matrix[,MODifieR_input$group_indici[[2]]])

  beta1= beta + value #20^sftR + value #user defined parameter for soft thresholding
  AdjMatC1<-sign(cor(dataset1, method = cor_method))*(cor(dataset1, method = cor_method))^2
  AdjMatC2<-sign(cor(dataset2 ,method = cor_method))*(cor(dataset2, method = cor_method))^2
  diag(AdjMatC1)<-0
  diag(AdjMatC2)<-0
  WGCNA::collectGarbage()

  dissTOMC1C2 <- TOMdist((abs(AdjMatC1-AdjMatC2)/2)^(beta1/2))
  WGCNA::collectGarbage()

  #Hierarchical clustering is performed using the Topological Overlap of the adjacency difference as input distance matrix
  geneTreeC1C2 <<- flashClust::flashClust(as.dist(dissTOMC1C2), method = cluster_method);

  # Plot the resulting clustering tree (dendrogram)
  #png(file="hierarchicalTree.png",height=1000,width=1000)
  #plot(geneTreeC1C2, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",labels = FALSE, hang = 0.04);
  #dev.off()

  #We now extract modules from the hierarchical tree. This is done using cutreeDynamic. Please refer to WGCNA package documentation for details
  dynamicModsHybridC1C2 <- dynamicTreeCut::cutreeDynamic(dendro = geneTreeC1C2,
                                                         distM = dissTOMC1C2,
                                                         method = cuttree_method,
                                                         cutHeight = cut_height,
                                                         deepSplit = deep_split,
                                                         pamRespectsDendro = pam_respects,
                                                         minClusterSize = min_cluster_size);
  
  #Every module is assigned a color. Note that GREY is reserved for genes which do not belong to any differentially coexpressed module
  dynamicColorsHybridC1C2 <<- labels2colors(dynamicModsHybridC1C2)

  #the next step merges clusters which are close (see WGCNA package documentation)
  colorh1C1C2 <<- mergeCloseModules(rbind(dataset1,dataset2),
                                     dynamicColorsHybridC1C2,
                                     cutHeight=cluster_merge_cut_height)$color
   #reassign better colors
  colorh1C1C2[which(colorh1C1C2 =="midnightblue")]<-"red"
  colorh1C1C2[which(colorh1C1C2 =="lightgreen")]<-"yellow"
  colorh1C1C2[which(colorh1C1C2 =="cyan")]<-"orange"
  colorh1C1C2[which(colorh1C1C2 =="lightcyan")]<-"green"
  # Plot the dendrogram and colors underneath
  #png(file="module_assignment.png",width=1000,height=1000)
  #plotDendroAndColors(geneTreeC1C2, colorh1C1C2, "Hybrid Tree Cut",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors cells")
  #dev.off()

  #We write each module to an individual file containing affymetrix probeset IDs
  modulesC1C2Merged <-extractModules(colorh1C1C2, dataset1, t(dataset1))
  
  annotation_table <- MODifieR_input$annotation_table
  
  module_list <- list()
  for (i in 1:length(modulesC1C2Merged)){
    module_list[[i]] <- construct_diffcoex_module(module = modulesC1C2Merged[[i]], annotation_table = annotation_table, 
                                                  module_color = names(modulesC1C2Merged)[[i]], settings = settings)
  }
  return(module_list)
}
