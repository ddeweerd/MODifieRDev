#' WGCNA
#'
#' A Coexpression based Algorithm to identify disease module
#'
#' @param dataset1 A dataframe of the microarray matrix of the Control samples
#' @param dataset2 A dataframe of the microarray matrix of the Patient samples
#' @param traits_file A dataframe of the associted traits/phenotypes  for the included samples

#' @export
wgcna <- function(MODifieR_input, pval_cutoff = 0.05, min_module_size = 30, deep_split = 2, pam_respects_dendro = F,
                  merge_cut_height = 0.1, numeric_labels = T, min_KME = 0,
                  save_toms = T, dataset_name = deparse(substitute(MODifieR_input))){
  # Retrieve settings
  default_args <- formals()
  user_args <- as.list(match.call(expand.dots = T)[-1])
  settings <- c(user_args, default_args[!names(default_args) %in% names(user_args)])
  settings$dataset_name <- NULL
  
  if (!is.null(dataset_name)){
    settings$MODifieR_input <- dataset_name
  }
  
  # The following setting is important, do not omit.
  options(stringsAsFactors = FALSE)
  allowWGCNAThreads()
  #Prepare data, traits
  Traits <- list()
  
  for (i in 1:length(MODifieR_input$group_indici)){
    trait <- data.frame("Label" = paste0(MODifieR_input$group_labels[[i]], 1:length(MODifieR_input$group_indici[[i]])),
                        "Treatment" = 1:length(MODifieR_input$group_indici[[i]]))
    rownames(trait) <- colnames(MODifieR_input$expression_matrix)[MODifieR_input$group_indici[[i]]]
    Traits[[i]] <- list("data" = trait)
  }
  setLabels <-  names(MODifieR_input$group_indici)
  #Prepare data, labels
  nSets = length(setLabels)
  
  #Prepare data, expression
  dataset1 <- list("data" = t(MODifieR_input$annotated_exprs_matrix[ ,MODifieR_input$group_indici[[1]]]))
  dataset2 <- list("data" = t(MODifieR_input$annotated_exprs_matrix[ ,MODifieR_input$group_indici[[2]]]))
  multiExpr <- list(dataset1,dataset2)
  
  names(multiExpr) <- c(paste(dataset_name, names(MODifieR_input$group_indici)[[1]], sep = "_"),
                        paste(dataset_name, names(MODifieR_input$group_indici)[[2]], sep = "_"))
  
  exprSize <- WGCNA::checkSets(multiExpr)
  
  gsg = WGCNA::goodSamplesGenesMS(multiExpr, verbose = 3)
  if (!gsg$allOK){
    # Print information about the removed genes:
    if (sum(!gsg$goodGenes) > 0)
      dynamicTreeCut::printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes],
                                                                collapse = ", ")))
    for (set in 1:exprSize$nSets)
    {
      if (sum(!gsg$goodSamples[[set]]))
        dynamicTreeCut::printFlush(paste("In set", setLabels[set], "removing samples",
                                         paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
      # Remove the offending genes and samples
      multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes]
    }
    # Update exprSize
    exprSize = WGCNA::checkSets(multiExpr)
  }
  
  WGCNA::collectGarbage()
  # Define data set dimensions
  nGenes = exprSize$nGenes
  nSamples = exprSize$nSamples
  
  # Choose a set of soft-thresholding powers
  powers = c(seq(1,15,by=1), seq(11,23, by=2))
  # Initialize a list to hold the results of scale-free analysis
  powerTables = vector(mode = "list", length = nSets)
  # Call the network topology analysis function for each set in turn
  power_modules <- NULL
  for (set in 1:nSets){
    powerTables[[set]] = list(data = WGCNA::pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                              verbose = 2)[[2]])
    
    x <- powerTables[[set]]$data[,1]
    y <- -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2]
    power_modules[set] <- min(x[which(y > max(y) * 0.95)])
  }
  print(names(multiExpr))
  net = WGCNA::blockwiseConsensusModules(
    multiExpr, power =  power_modules, minModuleSize = min_module_size,
    deepSplit = deep_split,
    pamRespectsDendro = pam_respects_dendro,
    mergeCutHeight = merge_cut_height,
    numericLabels = numeric_labels,
    minKMEtoStay = min_KME,
    saveTOMs = save_toms,
    individualTOMFileNames = paste0(dataset_name, "_individualTOM-Set%s-Block%b.RData"),
    verbose = 5)
  
  consMEs = net$multiMEs
  moduleLabels = net$colors
  # Convert the numeric labels to color labels
  moduleColors = WGCNA::labels2colors(moduleLabels)
  consTree = net$dendrograms[[1]];
  rm(net)
  
  
  # Set up variables to contain the module-trait correlations
  moduleTraitCor = list()
  moduleTraitPvalue = list()
  # Calculate the correlations
  for (set in 1:nSets)
  {
    moduleTraitCor[[set]] = WGCNA::cor(consMEs[[set]]$data, Traits[[set]]$data$Treatment, use = "p")
    moduleTraitPvalue[[set]] = WGCNA::corPvalueFisher(moduleTraitCor[[set]], exprSize$nSamples[set])
  }
  
  # Convert numerical lables to colors for labeling of modules in the plot
  MEColors = WGCNA::labels2colors(as.numeric(substring(names(consMEs[[1]]$data), 3)))
  MEColorNames = paste("ME", MEColors, sep="")
  
  # Initialize matrices to hold the consensus correlation and p-value
  consensusCor = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]))
  consensusPvalue = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]))
  
  # Find consensus negative correlations
  negative = moduleTraitCor[[1]] < 0 & moduleTraitCor[[2]] < 0
  consensusCor[negative] = pmax(moduleTraitCor[[1]][negative], moduleTraitCor[[2]][negative])
  consensusPvalue[negative] = pmax(moduleTraitPvalue[[1]][negative], moduleTraitPvalue[[2]][negative])
  # Find consensus positive correlations
  positive = moduleTraitCor[[1]] > 0 & moduleTraitCor[[2]] > 0
  consensusCor[positive] = pmin(moduleTraitCor[[1]][positive], moduleTraitCor[[2]][positive])
  consensusPvalue[positive] = pmax(moduleTraitPvalue[[1]][positive], moduleTraitPvalue[[2]][positive])
  
  row.names(consensusCor) <- MEColorNames
  new_consensusCor <- na.omit(consensusCor)

  
  selectiv_consensusCor <- data.frame(new_consensusCor[abs(new_consensusCor[,1]) > pval_cutoff,])
  
  correlation_colors <- selectiv_consensusCor
  
  consensusCor_up <- data.frame(new_consensusCor[(new_consensusCor[,1]>0),])
  consensusCor_down <- data.frame(new_consensusCor[(new_consensusCor[,1]<0),])
  colnames(consensusCor_up) <- NULL
  colnames(consensusCor_down) <- NULL
  colnames(selectiv_consensusCor) <- NULL
  
  probes <- colnames(multiExpr[[1]]$data)
  
  consMEs.unord = WGCNA::multiSetMEs(multiExpr, universalColors = moduleLabels, excludeGrey = TRUE)
  GS = list()
  kME = list()
  for (set in 1:nSets)
  {
    GS[[set]] = WGCNA::corAndPvalue(multiExpr[[set]]$data, Traits[[set]]$data$Treatment);
    kME[[set]] = WGCNA::corAndPvalue(multiExpr[[set]]$data, consMEs.unord[[set]]$data);
  }
  
  GS.metaZ = (GS[[1]]$Z + GS[[2]]$Z)/sqrt(2);
  kME.metaZ = (kME[[1]]$Z + kME[[2]]$Z)/sqrt(2);
  GS.metaP = 2*pnorm(abs(GS.metaZ), lower.tail = FALSE);
  kME.metaP = 2*pnorm(abs(kME.metaZ), lower.tail = FALSE);
  
  GSmat = rbind(GS[[1]]$cor, GS[[2]]$cor, GS[[1]]$p, GS[[2]]$p, GS.metaZ, GS.metaP);
  #nTraits = checkSets(Traits)$nGenes
  nTraits = 1
  traitNames = colnames(Traits[[1]]$data)
  dim(GSmat) = c(nGenes, 6*nTraits)
  
  rownames(GSmat) = probes
  colnames(GSmat) = WGCNA::spaste(
    c("GS.set1.", "GS.set2.", "p.GS.set1.", "p.GS.set2.", "Z.GS.meta.", "p.GS.meta"),
    rep(traitNames[2], rep(6, nTraits)))
  # Same code for kME:
  kMEmat = rbind(kME[[1]]$cor, kME[[2]]$cor, kME[[1]]$p, kME[[2]]$p, kME.metaZ, kME.metaP);
  MEnames = colnames(consMEs.unord[[1]]$data);
  nMEs = WGCNA::checkSets(consMEs.unord)$nGenes
  dim(kMEmat) = c(nGenes, 6*nMEs)
  rownames(kMEmat) = probes
  colnames(kMEmat) = WGCNA::spaste(
    c("kME.set1.", "kME.set2.", "p.kME.set1.", "p.kME.set2.", "Z.kME.meta.", "p.kME.meta"),
    rep(MEnames, rep(6, nMEs)))
  
  annotation_table <- MODifieR_input$annotation_table
  
  matched_probes <- unique(annotation_table$PROBEID) %in% probes
  
  info <- data.frame(Probe = probes, 
               ModuleLabel = moduleLabels,
               ModuleColor = WGCNA::labels2colors(moduleLabels),
               GSmat,
               kMEmat, stringsAsFactors = F)
  
  rownames(selectiv_consensusCor) <- gsub("ME", "", rownames(selectiv_consensusCor))
  module_genes <-info[info$ModuleColor %in% rownames(selectiv_consensusCor),]

  module_genes <- rownames(module_genes)
  module_genes <- module_genes[!is.na(module_genes)]
  
  new_wgcna_module <- list("module_genes" =  module_genes,
                           "probe_info_table" = info[, 1:3],
                           "module_correlation" = correlation_colors,
                           "settings" = settings)
  
  class(new_wgcna_module) <- c("MODifieR_module", "WGCNA")
  
  return (new_wgcna_module)
}
