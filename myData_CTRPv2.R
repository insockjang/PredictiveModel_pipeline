myData_CTRPv2 <- function(data.type = c("E","E_rnaseq","EL"),drug.type = c("ActArea","IC50")){
  
  library(predictiveModeling)
  library(synapseClient)
  synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
  
  #### Load CCLE Molecular Feature Data from Synapse ####    
  load("/gpfs/archive/RED/isjang/Project_CC90009/CCLE/eSet_expr.rdata")
  
  id_lineageLayer <- "syn1744846"  
  layer_lineage <- loadEntity(id_lineageLayer)
  eSet_lineage <- layer_lineage@objects$eSet_lineage
  
  myActArea<-function(){
    load("/gpfs/archive/RED/isjang/Project_CC90009/CTRP/CTRPv2.Rdata")
    adf_drug <- CTRPv2
    return(adf_drug)
  }
  
  # myIC50<-function(){
  #   load("/gpfs/archive/RED/isjang/Project_CC90009/PRISM/Celgene_IC50_decoded_jis.rdata")
  #   adf_drug <- PRISM.ic50
  #   return(adf_drug)
  # }
  
  drug.fun <- match.arg(drug.type)
  switch(drug.fun, 
         ActArea = (myfun1 = myActArea),
         IC50 = (myfun1 = myIC50))     
  adf_drug<-myfun1()
  
  myE <- function(){
    featureData <- createAggregateFeatureDataSet(list(expr = eSet_expr))    
    featureData_filtered <- filterNasFromMatrix(featureData, filterBy = "rows")
    dataSets <- createFeatureAndResponseDataList(t(featureData_filtered),t(adf_drug))
    return(dataSets)
  }
  myEranseq <- function(){
    featureData <- createAggregateFeatureDataSet(list(expr = eSet_expr))    
    featureData_filtered <- filterNasFromMatrix(featureData, filterBy = "rows")
    dataSets <- createFeatureAndResponseDataList(t(featureData_filtered),t(adf_drug))
    return(dataSets)
  }
  myEL <- function(){
    featureData <- createAggregateFeatureDataSet(list(expr = eSet_expr,line = eSet_lineage))    
    featureData_filtered <- filterNasFromMatrix(featureData, filterBy = "rows")
    dataSets <- createFeatureAndResponseDataList(t(featureData_filtered),t(adf_drug))
    return(dataSets)
  }
  
  data.fun <- match.arg(data.type)
  switch(data.fun, 
         E = (myfun = myE),
         EL = (myfun = myEL),
         E_rnaseq = (myfun = myErnaseq))     
  
  dataSets<-myfun()
  return(dataSets)
}
