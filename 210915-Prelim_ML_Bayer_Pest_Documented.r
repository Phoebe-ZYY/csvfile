############################################################################

### ToxCast/Tox21 active chem-endpoint associations 
###
### Author: Marissa Kosnik

############################################################################

rm(list = ls())

library(data.table)
library(prodlim)
library(caret)
library(randomForest)

#Test datasets

setwd("C:/me/map/cropmap")
#Existed: complete info including chem use for crop-country combos in EU 
exist.dat<- fread("combi_cntr_pest_eu_existed.csv")
head(exist.dat)
summary(exist.dat)
table(exist.dat$name_crop)
table(exist.dat$variable)
table(exist.dat$name_cntr)
table(exist.dat$name_ai)
table(exist.dat$layer_climate)
table(exist.dat$regu_eu)
table(exist.dat$regu_cntr)

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

summary_wdhaother

exist.dat.chem <- exist.dat %>%
  group_by(casrn_chemical) %>%
  mutate(count = n()) %>%
  filter(count >60) #%>%
  #filter(count < 10)

exist.dat.summary <- summary_wdhaother %>%
  filter(casrn_chemical %in% chem.list)

exist.dat.summary.country <- exist.dat.summary %>%
  group_by(name_country, name_region) %>%
  summarize(summass = sum(summass)) %>%
  left_join(country_continent, by = "name_country")
 write.csv(exist.dat.summary.country, "20220617chemicalwithfewscenarios_countrymass.csv")

summary.wdhaother.country %>%
  ggplot(aes(x = name_region.y, y = summass)) + 
  geom_bar(stat = "sum")

summary.wdhaother.country <- summary_wdhaother %>%
  group_by(name_country, name_region) %>%
  summarize(summass = sum(summass)) %>%
  left_join(country_continent, by = "name_country")
write.csv(summary.wdhaother.country, "20220617chemical_countrymass.csv")

 country_continent <- read_xlsx("C:/me/Database/harmonization/harmonization_country.xlsx")
# 
 sum(exist.dat.summary$summass)/sum(summary_wdhaother$summass)
# sum(exist.dat.chem$)
chem.list<- unique(exist.dat.chem$casrn_chemical)
unique(summary_wdhaother$casrn_chemical)

#Prep set of general scenarios that match chemical-specific scenarios
general.scenario<- exist.dat[, -c("casrn_chemical")]
general.scenario<- general.scenario[!duplicated(general.scenario)] #Lots of replicates in scenarios, remove

#Prep dataset to fill in with ML values
dat.out1<- matrix(data="X", nrow=length(chem.list), ncol= 12, dimnames=list(NULL, c("CASRN", "nNo", "nYes", "No_ActualNo", "Yes_ActualNo",  
                                                                                  "No_ActualYes", "Yes_ActualYes", "No_Accuracy", "Yes_Accuracy", "CPEfficiency", "CPErrorRate", "CPValidity")))
dat.out2<- matrix(data="X", nrow=length(chem.list), ncol= 12, dimnames=list(NULL, c("CASRN", "nNo", "nYes", "No_ActualNo", "Yes_ActualNo",  
                                                                                    "No_ActualYes", "Yes_ActualYes", "No_Accuracy", "Yes_Accuracy", "CPEfficiency", "CPErrorRate", "CPValidity")))
dat.out<- matrix(data="X", nrow=length(chem.list), ncol= 8, dimnames=list(NULL, c("CASRN", "nNo", "nYes", "No_Accuracy", "Yes_Accuracy", "CPEfficiency", "CPErrorRate", "CPValidity")))
#Building a separate model for each chem - loop through chem list and test accuracy
#for(i in c(1:length(chem.list))){
i <- 1
for(i in c(1:length(chem.list))){
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
  

  #Create a training set - this is what will be used to build the model (random 70% of machine data rows)
  seed <- 123
  set.seed(seed)
  train.dat<- createDataPartition(machine.dat$ChemPresent, p=0.7, list=FALSE)
  
  #The remaining 30% of rows are the test set - these will determine how good the model is
  test.dat<- machine.dat[-train.dat,] 
  
  train.dat<- machine.dat[train.dat,] 
  
  #Do random forest - just using base settings for now
  # rf <- randomForest(
  #   ChemPresent ~ .,
  #   data=train.dat
  # )
  
  pValues = ICPClassification(train.dat, test.dat, method = "rf",nrTrees = 100)
  pValues.df <- as.data.frame(pValues)
  testLabels = test.dat[,1]
  CPEfficiency <- CPEfficiency(pValues, testLabels)
  CPErrorRate <- CPErrorRate(pValues, testLabels)
  CPValidity <- CPValidity(pValues, testLabels)
  
  sig <- 0.15
  
  pValues.df.pred <- pValues.df %>%
    #mutate(chemPred = names(.)[max.col(.)])
    mutate(chemPred = case_when(V1 >= sig & V2 >= sig ~ "both p >= sig",
                                V1 < sig & V2 < sig ~ "both p < gis",
                                V1 >= sig & V2 < V1 ~ "V1",
                                V2 >= sig & V2 > V1 ~ "V2"))
  chemPred <- pValues.df.pred$chemPred
  #test.dat.p <- cbind(test.dat, pValues.df.pred)
  #test.dat.p$chemPred <-gsub("V","",as.character(test.dat.p$chemPred))
  #test.dat.p$ChemPresent <- as.character(test.dat.p$ChemPresent)
  #test.dat.p <- as.data.table(test.dat.p)
  #model <- fitModel(train.dat, method = "rf", nrTrees = 100)
  #Use the model to predict the test set
  # test <- test.dat %>%
  #   select(-ChemPresent)
  # pred <- predict(model, test)
  test.dat$ChemPresent <- as.character(test.dat$ChemPresent)
  #How well was the test set predicted? develop confusion matrix
  cm <- as.data.table(table(test.dat$ChemPresent, chemPred))
  
  #Put together data values to track accuracy of models
  val.out<- c(chem,
              nrow(absent.dat), #Total Nos in dataset
              nrow(present.dat), #Total yes in dataset
              cm[V2 =="1" & chemPred == "V1"]$N, #No actual no
              cm[V2 =="1" & chemPred == "V2"]$N, #Yes actual no
              cm[V2 =="2" & chemPred == "V1"]$N, #No actual yes
              cm[V2 =="2" & chemPred == "V2"]$N, #Yes actual yes
              cm[V2 =="1" & chemPred == "V1"]$N/length(which(test.dat$ChemPresent == "1")), #No accuracy
              cm[V2 =="2" & chemPred == "V2"]$N/length(which(test.dat$ChemPresent == "2")), #Yes accuracy
              CPEfficiency,
              CPErrorRate,
              CPValidity
              ) 
  # CPEfficiency <- round(CPEfficiency, 3)
  # CPErrorRate <- round(CPErrorRate, 3)
  # CPValidity <- round(CPValidity, 3)
  # 
  # jpeg(file = paste0("220323ML_", chem, ".jpeg"))
  # CPCalibrationPlot(pValues, test.dat, "blue")
  # title(main = paste0(chem," ", CPEfficiency," ", CPErrorRate, " ", CPValidity), cex = 0.1)
  # dev.off()
  
  dat.out1[i,]<- val.out
  print(i)
}


write.csv(dat.out1, "20220608dat.outwithoutcntrwithcropscorewithclimatemain20to60sig0.15.csv")

dat.out1 <- as.data.frame(dat.out1)
dat.out2 <- as.data.frame(dat.out2)
dat.out12 <- dat.out1 %>%
  full_join(dat.out2, by = c("CASRN", "nYes", "nNo"))

#x represents sig 0.15, y represents sig 0.05
dat.out12 <- dat.out12 %>%
  mutate(yes_diff = (as.numeric(Yes_Accuracy.x) - as.numeric(Yes_Accuracy.y))/as.numeric(Yes_Accuracy.y),
         no_diff = (as.numeric(No_Accuracy.x) - as.numeric(No_Accuracy.y))/as.numeric(No_Accuracy.y))
dat.out12 %>%
  ggplot(aes(x = yes_diff)) +
  geom_histogram()
write.csv(dat.out12, "20220608dat.out.sigfull2.csv")


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

dat.out.sig15 <- read.csv("C:/Users/yuyzh/Documents/20220606dat.outwithoutcntrwithcropscorewithclimatemain10to40.csv")
dat.out.sig5 <- read.csv("C:/Users/yuyzh/Documents/20220606dat.outwithoutcntrwithcropscorewithclimatemain10to40sig0.05.csv")
dat.out.sigfull <- dat.out.sig5 %>%
  full_join(dat.out.sig15, by = "CASRN")
write.csv(dat.out.sigfull, "20220608dat.out.sigfull.csv")

pred_df_108_62_3 <- read.csv("C:/me/data gaps/20220616pred_dffull108-62-3.csv")

pred_108_62_3 <- pred_df_108_62_3$pred
orig_108_62_3 <- pred_df_108_62_3$wd
rss <- sum((pred_108_62_3 - orig_108_62_3) ^ 2)  ## residual sum of squares
tss <- sum((orig_108_62_3 - mean(orig_108_62_3)) ^ 2)  ## total sum of squares
rsq <- 1 - rss/tss
rsq(df)

cor(orig_108_62_3, pred_108_62_3) ^ 2
