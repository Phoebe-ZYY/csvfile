############################################################################

### Application method machine learning
###
### Author: Marissa Kosnik

############################################################################

rm(list = ls())

library(data.table)
library(prodlim)
library(caret)
library(randomForest)
library(tidyverse)
library(grid)
library(gridExtra)
library(xgboost)
#library(ranger)

Sys.setenv(LANG = "en")

setwd("C:/Users/markos/Documents/Coding_Help/Bayer_Data_Gaps")
app.dat<- fread("C:/me/data gaps/application_input.csv")
head(app.dat)
str(app.dat) #See data types/levels

#Are there NAs in any predictors?
test<- app.dat[is.na(HDI_2015) == FALSE & is.na(GDP_2015) == FALSE & is.na(landevap_2017) == FALSE & 
                  is.na(tmin_all_mean) == FALSE & is.na(tmax_all_mean) == FALSE & is.na(precip_all_mean) == FALSE & 
                  is.na(ocs) == FALSE & is.na(soc_15_proj) == FALSE & is.na(soc_5_proj) == FALSE]

nrow(test)/nrow(app.dat) #work with this 83% for now
app.dat<- test #Remove rows with NA predictors

summary(app.dat)

#Remove duplicate values
length(unique(app.dat$name_ai.x)) == length(unique(app.dat$casrn_chemical)) #Only use 1 of these columns

app.dat<- app.dat[, -c("V1", "name_ai.x")]

#Check on different input data present
table(app.dat$name_refapp)

#Make sure no duplicates
app.dat<- app.dat[!duplicated(app.dat)]

#Fix input data types - factor for category, numeric for continuous
app.dat<- transform(
  app.dat,
  name_country=as.factor(as.character(name_country)),
  name_crop_earthstat=as.factor(as.character(name_crop_earthstat)),
  name_cropgroup.y=as.factor(as.character(name_cropgroup.y)),
  casrn_chemical=as.factor(as.character(casrn_chemical)),
  class1=as.factor(as.character(class1)),
  name_indication=as.factor(as.character(name_indication)),
  name_refapp=as.factor(as.character(name_refapp)),
  layer_climate=as.factor(as.character(layer_climate)),
  HDI_2015=as.numeric(as.character(HDI_2015)),
  GDP_2015=as.numeric(as.character(GDP_2015)),
  landevap_2017=as.numeric(as.character(landevap_2017)),
  gpw_2015=as.numeric(as.character(gpw_2015)),
  tmin_all_mean=as.numeric(as.character(tmin_all_mean)),
  tmax_all_mean=as.numeric(as.character(tmax_all_mean)),
  precip_all_mean=as.numeric(as.character(precip_all_mean)),
  ocs=as.numeric(as.character(ocs)),
  soc_15_proj=as.numeric(as.character(soc_15_proj)),
  soc_5_proj=as.numeric(as.character(soc_5_proj))
)


############################################################################

#Select columns with few enough categories
app.in<- app.dat[, list(name_country, name_cropgroup.y, name_indication, name_refapp, layer_climate,
                          HDI_2015, GDP_2015, landevap_2017, gpw_2015, tmin_all_mean, tmax_all_mean, precip_all_mean,
                          ocs, soc_15_proj, soc_5_proj)]
app.red<- app.in[sample(nrow(app.in), (nrow(app.in)/2))] #Downsample data

###Function to assess the accuracy of different ML methods
ml.bayer.test<- function(in.dat, test.var){
  
  test.var<- as.formula(paste(test.var, "~ ."))
  
  control<- trainControl(method="repeatedcv", number=5, repeats=2) #This will be the control for each model run
  
  fit.lda<- train(test.var, data=in.dat, method="lda", metric="Accuracy", trControl=control)
  
  fit.cart<- train(test.var, data=in.dat, method="rpart", metric="Accuracy", trControl=control)
  
  #fit.knn<- train(test.var, data=in.dat, method="knn", metric="Accuracy", trControl=control)
  
  #fit.svm<- train(test.var, data=in.dat, method="svmRadial", metric="Accuracy", trControl=control)
  
  fit.rf<- train(test.var, data=in.dat, method="ranger", metric="Accuracy", trControl=control)
  
  #fit.nn<- train(test.var, data=in.dat, method="nnet", metric="Accuracy", trControl=control)
  
  #fit.gbm<- train(test.var, data=in.dat, method="gbm", metric="Accuracy", trControl=control)
  
  fit.xgb<- train(test.var, data=in.dat, method="xgbTree", metric="Accuracy", trControl=control)
  
  # summarize accuracy of models
  dat.out <- resamples(list(lda=fit.lda, cart=fit.cart, rf=fit.rf, xgb=fit.xgb))
  dat.out
}

app.test<- ml.bayer.test(app.red, "name_refapp") 
ysummary(app.test) #Good accuracy across the board, best with xgboost - moving forward with this data

###Now need to train and test model

#2.3 in walkthrough
train.dat<- createDataPartition(app.red$name_refapp, p=0.7, list=FALSE) #Select random rows to be the training data
test.dat<- app.red[!train.dat] #Select other rows to be the test data
train.dat<- app.red[train.dat]

#Using xgboost because it had the best accuracy (this plus the function above = step 5 in walkthrough)

control<- trainControl(method="repeatedcv", number=5, repeats=2) #This will be the control for each model run
fit.app.model<- train(name_refapp~., data=train.dat, method="xgbTree", metric="Accuracy", trControl=control)

#Make predictions on test set - Step 6 in walkthrough
predict.app<- predict(fit.app.model, test.dat)

#Look at confusion matrix
cm<- confusionMatrix(predict.app, test.dat$name_refapp) #This compares predictions to actual values
cm.out<- cm$table
cm.vals<- cm$byClass

fwrite(as.data.table(cm.out, keep.rownames = TRUE), file="211215-CM_AppMethod.csv")
fwrite(as.data.table(cm.vals, keep.rownames = TRUE), file="211215-Accuracy_AppMethod.csv")
#save.image("Train-Test_Application_Method.RData")
