Drug.type = c("ActArea")
Data.type = c("E","EL")
Model.type = c("Lasso","ENet")

k3 <- c(1)
source("myBSModel1.R")
source("myData_PRISM.R")

for(model.Type in Model.type){  
  for(drug.Type in Drug.type){
    for(kk in k3){
      for(data.Type in Data.type){
        filename = paste("~/newPredictiveModel/PRISM/",drug.Type,"/",data.Type,"/bs_",model.Type,"/bsDrug_",kk,".Rdata",sep="")
        if(!file.exists(filename)){
          resultsScale <- myBSModel1(kk,data.set = "PRISM",data.type=data.Type, drug.type = drug.Type, model.type = model.Type, numBS= 100, numCore = 4)    
          save(resultsScale,file = filename)        
        }  
      }  
    }
  }
}


