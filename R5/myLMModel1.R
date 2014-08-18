myLMModel1 <- setRefClass(Class = "myLMModel1",
                            contains="PredictiveModel",
                            fields="model",
                            methods = list(
                              initialize = function(...){
                                return(.self)
                              },
                              
                              rawModel = function(){
                                return(.self$model)
                              },
                              
                              customTrain = function(featureData, responseData,...){
                                
                                  .self$model <- lm(responseData ~ featureData,...)                                   
                                
                              },
                              
                              customPredict = function(featureData){
                                predictedResponse <- predict(.self$model, featureData)
                                return(predictedResponse)
                              },
                              
                              getCoefficients = function(){
                                return(coef(.self$model))
                              }
                              
                            )
                            
)
