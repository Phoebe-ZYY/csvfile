############################################################################

### chemical_imputation_seperate_extrapolation
###
### Author: Yuyue Zhang

############################################################################

rm(list = ls())

library(data.table)
library(prodlim)
library(caret)
library(randomForest)
library(dplyr)

#Test datasets

exist.dat <- fread("C:/me/data gaps/20220522country_input_full.csv")
exist.dat$climate_main <- substr(exist.dat$layer_climate, 1, 1)
exist.dat$climate_sub <- substr(exist.dat$layer_climate, 2, 2)


#Existing countries with missing crops 
nocrop.dat<- fread("combi_cntr_pest_eu_nocrop.csv")
head(nocrop.dat)
summary(nocrop.dat)
table(nocrop.dat$name_crop)
table(nocrop.dat$variable)
table(nocrop.dat$name_cntr)
table(nocrop.dat$name_ai)
table(nocrop.dat$layer_climate)
table(nocrop.dat$regu_eu)
table(nocrop.dat$regu_cntr)
length(which(!nocrop.dat$name_crop %in% exist.dat$name_crop))
length(which(!unique(nocrop.dat$name_crop) %in% exist.dat$name_crop))

#Test countries with no chem crop info - adding banned chems
#norecord.dat<- fread("combi_cntr_pest_eu_norecord.csv")

#First goal: Predict for countries with some data
#First: just see how train/test accuracy is with existing data

#Check data types again
summary(exist.dat)


#Remove unusable values
length(which(is.na(exist.dat$regu_eu)==TRUE))
length(which(exist.dat$name_cntr != exist.dat$name_country)) #Only need 1 of these columns
length(unique(exist.dat$variable)) == length(unique(exist.dat$name_crop)) #Only use 1 of these columns
length(unique(exist.dat$casrn_chemical)) == length(unique(exist.dat$name_ai)) #Only use 1 of these columns

exist.dat<- exist.dat[, -c("V1", "summass", "sumarea", "mindose", "name_chemclass", "wd", "name_country")]
#exist.dat<- exist.dat[, -c("count")]

exist.dat <- exist.dat[is.na(HDI) == FALSE & is.na(landevap) == FALSE]

chem.list<- unique(exist.dat$casrn_chemical)

#Prep set of general scenarios that match chemical-specific scenarios
general.scenario<- exist.dat[, -c("casrn_chemical")]
general.scenario<- general.scenario[!duplicated(general.scenario)] #Lots of replicates in scenarios, remove

#Prep dataset to fill in with ML values
dat.out1<- matrix(data="X", nrow=length(chem.list), ncol= 9, dimnames=list(NULL, c("CASRN", "nNo", "nYes", "No_ActualNo", "Yes_ActualNo",  
                                                                                    "No_ActualYes", "Yes_ActualYes", "No_Accuracy", "Yes_Accuracy")))

#for(i in c(1:length(chem.list))){
for(i in c(1160:length(chem.list))){
  chem<- chem.list[i]
  present.dat<- exist.dat[casrn_chemical %in% chem] #ID scenarios that use that chem
  present.dat<- present.dat[, -c("casrn_chemical")]
  present.dat<- present.dat[!duplicated(present.dat)]
  x<- row.match(present.dat, general.scenario) #ID scenarios that do not use that chem
  absent.dat<- general.scenario[!na.omit(x)]
  
  #Do any of the present rows match absent rows? If yes, there's an error
  test<- row.match(present.dat, absent.dat)
  if(length(which(is.na(test) == FALSE)) != 0) stop ('Multi-value match')
  
  #Does a scenario use a chem?
  present.dat$ChemPresent<- "Yes"
  absent.dat$ChemPresent<- "No"
  
  machine.dat<- as.data.table(rbind(present.dat, absent.dat))
  write.csv(machine.dat,"machine.dat.csv")
  machine.dat <- read.csv("machine.dat.csv")
  # #Set machine dat types to match necessary input for ML 
  # machine.dat<- transform(
  #   machine.dat,
  #   ChemPresent=as.factor(ChemPresent),
  #   layer_climate=as.factor(layer_climate),
  #   #fieldsize_proj=as.integer(fieldsize_proj),
  #   HDI=as.numeric(HDI),
  #   #GDP_2015=as.numeric(GDP_2015),
  #   landevap=as.numeric(landevap),
  #   #gpw_2015=as.numeric(gpw_2015),
  #   #tmin_all_mean=as.numeric(tmin_all_mean),
  #   #tmax_all_mean=as.numeric(tmax_all_mean),
  #   #precip_all_mean=as.numeric(precip_all_mean),
  #   #ocs=as.numeric(ocs),
  #   #soc_15_proj=as.numeric(soc_15_proj),
  #   #soc_5_proj=as.numeric(soc_5_proj),
  #   name_country=as.factor(name_country),
  #   name_crop_earthstat=as.factor(name_crop_earthstat),
  #   name_region = as.factor(name_region),
  #   name_cropgroup = as.factor(name_cropgroup)
  # )
  machine.dat <- machine.dat %>%
    select(ChemPresent, name_crop_earthstat, name_region, name_cropgroup,HDI, landevap, layer_climate,climate_main, climate_sub,score)
  
  ## make sure first column is always the label and class labels are always 1, 2, ...
  machine.dat[, 1] = as.factor(machine.dat[, 1])
  machine.dat[, 1] = as.numeric(unlist(machine.dat[, 1]))
  machine.dat <- machine.dat %>%
    mutate(ChemPresent = case_when(ChemPresent == 1 ~ 0,
                                   ChemPresent == 2 ~ 1)) %>%
    select(ChemPresent, name_crop_earthstat, name_region, name_cropgroup,HDI, landevap, layer_climate,climate_main, climate_sub,score)
  machine.dat[, 1] = as.factor(machine.dat[, 1])
  
  
  #Create a training set - this is what will be used to build the model (random 70% of machine data rows)
  train.dat<- createDataPartition(machine.dat$ChemPresent, p=0.7, list=FALSE)
  
  #The remaining 30% of rows are the test set - these will determine how good the model is
  test.dat<- machine.dat[-train.dat,] 
  
  train.dat<- machine.dat[train.dat,] 
  
  #Do random forest - just using base settings for now
  # rf <- randomForest(
  #   ChemPresent ~ .,
  #   data=train.dat
  # )
  
  mymodel<-glm(ChemPresent~name_crop_earthstat + name_region + HDI + landevap + layer_climate + score, 
               data = train.dat, 
               family = 'binomial')
  
  ChemPred <- predict(mymodel, machine.dat,
                      type = 'response')
  Pred <- ifelse(ChemPred > 0.5, 1, 0)
  
  #test.dat.p <- cbind(test.dat, pValues.df.pred)
  #test.dat.p$chemPred <-gsub("V","",as.character(test.dat.p$chemPred))
  #test.dat.p$ChemPresent <- as.character(test.dat.p$ChemPresent)
  #test.dat.p <- as.data.table(test.dat.p)
  #model <- fitModel(train.dat, method = "rf", nrTrees = 100)
  #Use the model to predict the test set
  # test <- test.dat %>%
  #   select(-ChemPresent)
  # pred <- predict(model, test)
  machine.dat$ChemPresent <- as.character(machine.dat$ChemPresent)
  #How well was the test set predicted? develop confusion matrix
  cm <- as.data.table(table(machine.dat$ChemPresent, Pred))
  
  #Put together data values to track accuracy of models
  val.out<- c(chem,
              nrow(absent.dat), #Total Nos in dataset
              nrow(present.dat), #Total yes in dataset
              cm[V2 =="0" & Pred == "0"]$N, #No actual no
              cm[V2 =="0" & Pred == "1"]$N, #Yes actual no
              cm[V2 =="1" & Pred == "0"]$N, #No actual yes
              cm[V2 =="1" & Pred == "1"]$N, #Yes actual yes
              cm[V2 =="0" & Pred == "0"]$N/length(which(machine.dat$ChemPresent == "0")), #No accuracy
              cm[V2 =="1" & Pred == "1"]$N/length(which(machine.dat$ChemPresent == "1")) #Yes accuracy
  )
  
  dat.out1[i,]<- val.out
  print(i)
}

summary(mymodel)

write.csv(dat.out1, "20220530dat.outextrawithoutcntrwithcropscorewithclimatemain.csv")

#compare extrapolation result & cp result
dat.out.cp <- read.csv("20220523dat.outwithoutcntrwithcropscorewithclimatemain.csv")
dat.out1 <- as.data.frame(dat.out1)
dat.out.compare <- dat.out.cp %>%
  full_join(dat.out1, by = c("CASRN"))
write.csv(dat.out.compare, "20220531dat.outcomparewithoutcntrwithcropscorewithclimatemain.csv")

#check for first chemical, chempresent 2 for yes, 1 for no
pValues.df <- pValues.df %>%
  mutate(chemPred = names(.)[max.col(.)])
test.dat.p <- cbind(test.dat, pValues.df)
test.dat.p$chemPred <-gsub("V","",as.character(test.dat.p$chemPred))



jpeg(file="saving_plot1.jpeg")
CPCalibrationPlot(pValues, test.dat, "blue")
title(main = "plot1", cex= 0.1)
dev.off()
p
model$confusion
model$confusion[1,1]
save.image(file=paste0("220323ML_", chem, ".jpg")) 
dat.out<- as.data.table(dat.out)
write.csv(dat.out, "220323dat.out1.csv")
dat.out$No_Accuracy<- as.numeric(as.character(dat.out$No_Accuracy))
dat.out$Yes_Accuracy<- as.numeric(as.character(dat.out$Yes_Accuracy))
dat.out$nNo<- as.numeric(as.character(dat.out$nNo))
dat.out$nYes<- as.numeric(as.character(dat.out$nYes))

range(dat.out$No_Accuracy); range(dat.out$Yes_Accuracy)
quantile(dat.out$No_Accuracy); quantile(dat.out$Yes_Accuracy)

#Data viz of results

plot(dat.out$nNo, dat.out$No_Accuracy)
plot(dat.out$nYes, dat.out$Yes_Accuracy)
plot(dat.out$No_Accuracy, dat.out$Yes_Accuracy)

setwd("C:/Users/markos/Documents/DataDownloads/combi_cntr_crop_pest")
#save.image(file="210915-Preliminary_ML_Analysis.RData")
save.image(file="210916-Run2_Preliminary_ML_Analysis.RData")
#load("210915-Preliminary_ML_Analysis.RData")
