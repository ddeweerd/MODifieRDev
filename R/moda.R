#' MODA
#'
#'A Coexpression based algorithm to identify disease module
#' 
#' @param datExpr1 A dataframe of the microarray matrix of the Control samples
#' @param datExpr2 A dataframe of the microarray matrix of the Control samples
#' @export
moda <- function(MODifieR_input,
                 CuttingCriterion = "Density",
                 specificTheta = 0.1, conservedTheta = 0.1,
                 result_folder, dataset_name = NULL){
  # Retrieve settings
  default_args <- formals()
  user_args <- as.list(match.call(expand.dots = T)[-1])
  settings <- c(user_args, default_args[!names(default_args) %in% names(user_args)])
  
  if (!is.null(dataset_name)){
    settings$MODifieR_input <- dataset_name
  }
  
  datExpr1 <- t(MODifieR_input$annotated_exprs_matrix[,MODifieR_input$group_indici[[1]]])
  datExpr2 <- t(MODifieR_input$annotated_exprs_matrix[,MODifieR_input$group_indici[[2]]])
  
  
  
  indicator1 = names(MODifieR_input$group_indici)[1]     # indicator for data profile 1
  indicator2 = names(MODifieR_input$group_indici)[2]  # indicator for data profile 2
  specificTheta = specificTheta #threshold to define condition specific modules
  conservedTheta = conservedTheta #threshold to define conserved modules
  
  dir.create(result_folder, showWarnings = F)
  
  intModules1 <- MODA::WeightedModulePartitionHierarchical(datExpr = datExpr1, foldername = result_folder, 
                                                           indicatename = indicator1, cutmethod = CuttingCriterion)
  
  intModules2 <- MODA::WeightedModulePartitionHierarchical(datExpr2,result_folder,
                                                           indicator2,CuttingCriterion)
  print("WMPD Indicators 1 done")
  
  MODA::CompareAllNets(ResultFolder = result_folder, intModules = intModules1, indicator =  indicator1, 
                       intconditionModules = intModules2, conditionNames = indicator2, specificTheta = specificTheta, 
                       conservedTheta = conservedTheta)
  setwd(paste(result_folder, indicator2, sep = "/"))
  sep_moduleIDs = read.table("sepcificModuleid.txt", sep = "\n");
  cons_moduleIDs = read.table("conservedModuleid.txt", sep = "\n");
  
  setwd(result_folder)
  modules_control_sep = list()
  modules_patient_sep = list()
  for (i in 1: length(t(sep_moduleIDs))){
    if(file.exists(paste("DenseModuleGene_control_",sep_moduleIDs[i,],".txt",  sep =""))){
      if (file.info(paste("DenseModuleGene_control_",sep_moduleIDs[i,],".txt",  sep =""))$size != 0){
        modules_control_sep[i] = read.table(paste("DenseModuleGene_control_",sep_moduleIDs[i,],".txt",  sep =""), sep = "\n" , );
      }}
    if(file.exists(paste("DenseModuleGene_patient_",sep_moduleIDs[i,],".txt",  sep =""))){
      if (file.info(paste("DenseModuleGene_patient_",sep_moduleIDs[i,],".txt",  sep =""))$size != 0){
        modules_patient_sep[i] = read.table(paste("DenseModuleGene_patient_",sep_moduleIDs[i,],".txt",  sep =""), sep = "\n");
      }}
  }
  print("Modules patients and controls set")
  modules_control_cons = list()
  modules_patient_cons = list()
  for (i in 1: length(t(cons_moduleIDs))){
    if(file.exists(paste("DenseModuleGene_control_",cons_moduleIDs[i,],".txt",  sep =""))){
      if (file.info(paste("DenseModuleGene_control_",cons_moduleIDs[i,],".txt",  sep =""))$size != 0){
        modules_control_cons[i] = read.table(paste("DenseModuleGene_control_",cons_moduleIDs[i,],".txt",  sep =""), sep = "\n" , );
      }}
    if(file.exists(paste("DenseModuleGene_patient_",cons_moduleIDs[i,],".txt",  sep =""))){
      if (file.info(paste("DenseModuleGene_patient_",cons_moduleIDs[i,],".txt",  sep =""))$size != 0){
        modules_patient_cons[i] = read.table(paste("DenseModuleGene_patient_",cons_moduleIDs[i,],".txt",  sep =""), sep = "\n");
      }}
  }
 
    ## Union of submodules
    sep_module_control = unique(as.character(unlist(modules_control_sep)))
    sep_module_patient = unique(as.character(unlist(modules_patient_sep)))
    
    
    overlap_sep = VennDiagram::calculate.overlap(list(sep_module_control, sep_module_patient)) # check overlap
    overlap_sep2 = (t(overlap_sep$a3))
    
    cons_module_control = unique(as.character(unlist(modules_control_cons)))
    cons_module_patient = unique(as.character(unlist(modules_patient_cons)))
    print("Union of submodules 1")
    
    overlap_cons = VennDiagram::calculate.overlap(list(cons_module_control, cons_module_patient)) # check overlap
    overlap_cons2 = (t(overlap_cons$a3))
    print("Union of submodules 2")
    total_overlap = unique(c(overlap_cons2, overlap_sep2))
  
  

  
  new_moda_module <- list("module_genes" =  total_overlap,
                          "settings" = settings)
  
  class(new_moda_module) <- c("MODifieR_module", "MODA")
  
  return (new_moda_module)
  
  
}
