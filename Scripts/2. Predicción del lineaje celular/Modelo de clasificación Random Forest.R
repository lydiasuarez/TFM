library("e1071")
library("caret")
library("randomForest")
library("parallel")
library("doParallel")
library("devtools")

all_dataframe = readRDS("./all_dataframe_Immune.rds")

# Particiones

trainIndex <- caret::createDataPartition(all_dataframe$BrainRegion,
                                         p=0.80,
                                         list = F)
data.train <- all_dataframe[trainIndex,]
data.test <- all_dataframe[-trainIndex,]


# Variables 

variable_salida = c("BrainRegion")
variables_entrada = setdiff(names(all_dataframe),variable_salida)


# Parallel

intervalStart <- Sys.time()
cluster <- makeCluster(20) 
registerDoParallel(cluster)


# Training

cat("Empieza el training")

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


system.time(
  zheisel.Immune.fit.rds<-train(data.train[variables_entrada],
                                data.train[[variable_salida]],
                                method='rf',
                                trControl = ctrl,
                                metric="Accuracy",
                                seeds=seeds))


saveRDS(zheisel.modelo.3br.fit, file = "./results_immune/zheisel.Neurons.fit.rds")
saveRDS(ctrl,file="./results_immune/ctrl.rds")

stopCluster(cluster)
registerDoSEQ()
