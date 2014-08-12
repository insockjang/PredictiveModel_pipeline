myModel1<-function(kk,
                  data.set = c("CCLE","Sanger"),
                  data.type = c("Mh","C","CMo","CMh","E","EMo","EMh","EC","ECMo","ECMh","MhL","CL","CMoL","CMhL","EL","EMoL","EMhL","ECL","ECMoL","ECMhL"), 
                  drug.type = c("ActArea","IC50","EC50"), 
                  model.type = c("ENet","Lasso","Ridge","RF","PCR","PLS","SVM"), 
                  nfolds = 5){
  require(predictiveModeling)
  require(synapseClient)
  synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
  source("~/PredictiveModel_pipeline/R5/crossValidatePredictiveModel1.R")
  
  myCCLE<-function(X,Y){
    dataSets<-myData_CCLE_new(X,Y)
    return(dataSets)
  }
  mySanger<-function(X,Y){
    dataSets<-myData_Sanger(X,Y)
    return(dataSets)
  }
  
  myENet<-function(X,Y){
    source("~/PredictiveModel_pipeline/R5/myEnetModel1.R")
    alphas =unique(createENetTuneGrid()[,1])    
    CV<-crossValidatePredictiveModel1(X, Y, model = myEnetModel1$new(), alpha = alphas, numFolds = nfolds)
    return(CV)
  }
  myLasso<-function(X,Y){
    source("~/PredictiveModel_pipeline/R5/myEnetModel1.R")        
    CV<-crossValidatePredictiveModel1(X, Y, model = myEnetModel1$new(), alpha = 1, numFolds = nfolds)
    return(CV)
  }
  myRidge<-function(X,Y){
    source("~/PredictiveModel_pipeline/R5/myEnetModel1.R")    
    
    library(corpcor)    
    GetRidgeGrid <- function(d, U, V, y, beta.max = 0.001, epsilon = 1e-6,
                             K = 100, upper.bound = TRUE) {
      r <- length(d)
      aux <- V %*% t(U) %*% matrix(y,ncol=1)
      if (upper.bound) {
        lambda.max <- max(d * abs(aux[1:r,1]))/beta.max
      }
      else {
        lambda.max <- max(d * abs(aux[1:r,1]))/beta.max - max(d^2)
      }
      exp(seq(log(lambda.max), log(epsilon * lambda.max), length.out = K))
    }
    
    
    fsvd <- fast.svd(X)
    
    lambda.grid <- GetRidgeGrid(d = fsvd$d, U = fsvd$u, V = fsvd$v, y=Y, K = 100)
        
    CV<-crossValidatePredictiveModel1(X, Y, model = myEnetModel1$new(), alpha = 0, numFolds = nfolds, lambda = lambda.grid)
    return(CV)
  }
  myRF<-function(X,Y){
    source("~/PredictiveModel_pipeline/R5/myRandomForestModel1.R")    
    CV<-crossValidatePredictiveModel1(X, Y, model = myRandomForestModel1$new(), ntree = 500)
    return(CV)
  }
  myPCR<-function(X,Y){
    require(pls)
    source("~/PredictiveModel_pipeline/R5/myPcrModel1.R")    
    CV<-crossValidatePredictiveModel1(X, Y, model = myPcrModel1$new(), ncomp=10)
    return(CV)
  }
  myPLS<-function(X,Y){
    require(pls)
    source("~/PredictiveModel_pipeline/R5/myPlsModel1.R")    
    CV<-crossValidatePredictiveModel1(X, Y, model = myPlsModel1$new(), ncomp=10)
    return(CV)
  }
  mySVM<-function(X,Y){    
    myTrControl=trainControl(method = "cv", number = 5, returnResamp = "all", verboseIter = TRUE)    
    CV<-crossValidatePredictiveModel1(X, Y, model = CaretModel$new(modelType = "svmLinear"), trControl = myTrControl, numFolds = nfolds)
    return(CV)
  }
  
  model.fun <- match.arg(model.type)
  switch(model.fun, 
         ENet = (myfun = myENet),
         Lasso = (myfun = myLasso),
         Ridge = (myfun = myRidge),
         RF = (myfun = myRF),
         SVM = (myfun = mySVM),
         PCR = (myfun = myPCR),
         PLS = (myfun = myPLS))
  
  set.fun <- match.arg(data.set)
  switch(set.fun, 
         CCLE = (myfun2 = myCCLE),
         Sanger = (myfun2 = mySanger))
         
  
  dataSet<-myfun2(data.type,drug.type)
  
  # data preprocessing for preselecting features
  filteredData<-filterPredictiveModelData(dataSet$featureData,dataSet$responseData[,kk,drop=FALSE])
  
  # filtered feature and response data
  filteredFeatureData  <- filteredData$featureData
  filteredFeatureData  <- t(unique(t(filteredFeatureData)))
  filteredResponseData <- filteredData$responseData
  
  ## scale these data    
  filteredFeatureDataScaled <- scale(filteredFeatureData)
  filteredResponseDataScaled <- scale(filteredResponseData)  
  
  resultsScale<-myfun(filteredFeatureDataScaled,filteredResponseDataScaled)
  
  return(resultsScale) 
}
