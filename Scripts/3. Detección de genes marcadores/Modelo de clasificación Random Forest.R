library("e1071")
library("caret")
library("randomForest")
library("parallel")
library("doParallel")
library("devtools")


brainregions = c("amygdala","cortex2","hippocampus","hypothalamus","midbrainventral","midbraindorsal","olfactory","spinalcord","striatumventral","thalamus")
celltypes = c("Astrocytes","Immune","Neurons","Oligos")

# Parallel

intervalStart <- Sys.time()
cluster <- makeCluster(20) 
registerDoParallel(cluster)


# Random Forest

for (region in brainregions) {
  for (cell in celltypes){
    all_dataframe = readRDS(paste0("./",region,"/",cell,"/",region,".",cell,".dataset.rf.rds"))
    
    cat ("Comienza el Random Forest para: ",region, "-", cell,"\n")
    trainIndex <- caret::createDataPartition(all_dataframe$brainregion,
                                             p=0.80,
                                             list = F)
    
    data.train <- all_dataframe[trainIndex,]
    data.test <- all_dataframe[-trainIndex,]
    
    
    # Variables 
    
    variable_salida = c("brainregion")
    variables_entrada = setdiff(names(all_dataframe),variable_salida)
    
    
    # Training
    
    cat("Empieza el training","\n")
    
    ctrl <- trainControl(method = "repeatedcv",
                         number = 3,
                         repeats= 1,
                         verboseIter = TRUE,
                         returnResamp = "final",
                         allowParallel = TRUE)
    
    set.seed(123456)
    seeds <- vector(mode = "list", length = 11)
    for(i in 1:5) seeds[[i]]<- sample.int(n=1000, 10)
    seeds[[11]]<-sample.int(1000, 1)
    
    cat("Empieza el modelo","\n")
    
    system.time(
      zheisel.train.fit.rds<-train(data.train[variables_entrada],
                                   data.train[[variable_salida]],
                                   method='rf',
                                   trControl = ctrl,
                                   metric="Accuracy",
                                   seeds=seeds))
    
    
    saveRDS(zheisel.train.fit.rds, file = paste0("./",region,"/",cell,"/",region,".",cell,".train.fit.rds"))
    saveRDS(data.test,file=paste0("./",region,"/",cell,"/",region,".",cell,".data.test.rds"))
  }
}


stopCluster(cluster)
registerDoSEQ()