myPlsModel1 <- setRefClass(Class = "myPlsModel1",
                           contains="PredictiveModel",
                           fields="model",
                           methods = list(
                             initialize = function(...){
                               return(.self)
                             },
                             
                             rawCaretModel = function(){
                               return(.self$model)
                             },
                             
                             customTrain = function(featureData, responseData, ncomp,...){
                               A<-list(response = responseData,obsMat = featureData)
                               #                               A1<-as.data.frame(A)
                               .self$model<- mvr(response ~ obsMat,ncomp,data = A,method = pls.options()$plsralg)                                       
                             },
                             
                             customPredict = function(featureData){
                               predictedResponse <- predict(.self$model, ncomp = .self$model$ncomp, newdata=featureData)
                               return(predictedResponse)
                             }
                           )
)
