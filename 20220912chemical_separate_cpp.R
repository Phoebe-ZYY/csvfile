############################################################################

### Chem imputation crop-country-chemical-climate
###
### Author: Yuyue Zhang

############################################################################

library(data.table)
#library(scrambler)
library(prodlim)
library(caret)
library(randomForest)
library(conformalClassification)
library(splitTools)
library(survival)
library(dplyr)
library(psych)
library(intsurv)
library(parsnip)
library(recipes)
library(rsample)
library(readr)
library(workflows)
library(npreg)
#library(ProbBayes)
library(tibble)
library(tidyr)
library(purrr)
#Test datasets

#Existed: complete info including chem use for crop-country combos
exist.dat<- fread("C:/me/data gaps/20221019country_input_full_all.csv")
#exist.dat<- fread("C:/me/data gaps/20220912country_input_full_all.csv")
exist.dat.nocrop <- fread("C:/me/data gaps/20221019country_imputation_combi_existcountrynocrop.csv")
#exist.dat.nocrop <- fread("C:/me/data gaps/20221006country_imputation_combi_existcountrynocrop.csv")
exist.dat <- exist.dat %>%
  select(-name_group_yy) %>%
  mutate(name_region.y = case_when(name_country == "CENTRAL AMERICA-CARIBBEAN" ~ "North America",
                                   TRUE ~ name_region.y)) %>%
  left_join(crop_id, by = c("name_crop_earthstat" = "CROPNAME")) %>%
  left_join(country_croparea, by = c("name_country", "name_crop_earthstat" = "variable")) %>%
  left_join(no.country.continent, by = c("name_country"))

# no.dat.missingfdi <- no.dat %>%
#   filter(is.na(GCF) | is.na(FDI_in)) %>%
#   distinct(name_country)

no.country.continent <- read_xlsx("C:/me/Database/harmonization/no.country.continent.xlsx", sheet = 1)
#no.dat <- fread("C:/me/data gaps/20220912country_input_n_ccc.csv")
no.dat <- fread("C:/me/data gaps/20221019country_input_n_all.csv")
#no.dat.country <- unique(no.dat$name_country)
#write.csv(no.dat.country, "no.dat.country.csv")
no.dat <- no.dat %>%
  filter(score != 0) %>%
  #filter(name_country %in% country_noprediction$name_country) %>%
  left_join(crop_id, by = c("variable" = "CROPNAME")) %>%
  mutate(name_group_yy = name_group_yy)%>%
  left_join(no.country.continent, by = c("name_country" = "name_country")) %>%
  left_join(country_croparea, by = c("name_country", "variable")) %>%
  mutate(name_region.y = name_region) %>%
  mutate(name_crop_earthstat = variable) %>%
  #select("name_group_yy", "name_crop_earthstat", "name_region.y","layer_climate","score", "HDI", "landevap", "croparea", "sumcroparea" ) %>%
  filter(!is.na(HDI)) %>%
  filter(!is.na(landevap))
#no.dat.country <- no.dat %>%
  #group_by(name_country, name_region.y) %>%
  #summarize(count = n())
no.dat$climate_main <- substr(no.dat$layer_climate, 1, 1)
no.dat$climate_sub <- substr(no.dat$layer_climate, 2, 2)
no.dat <- no.dat[!duplicated(no.dat)]
exist.dat.nocrop <- exist.dat.nocrop %>%
  filter(score != 0) %>%
  select("layer_climate", "name_country", "score", "variable", "HDI", "landevap") %>%
  left_join(crop_id, by = c("variable" = "CROPNAME")) %>%
  mutate(name_group_yy = name_group_yy)%>%
  left_join(no.country.continent, by = c("name_country" = "name_country")) %>%
  left_join(country_croparea, by = c("name_country", "variable")) %>%
  mutate(name_region.y = name_region) %>%
  mutate(name_crop_earthstat = variable) %>%
  select("name_group_yy", "name_crop_earthstat", "name_region.y","layer_climate","score", "HDI", "landevap", "croparea", "sumcroparea" ) %>%
  filter(!is.na(HDI)) %>%
  filter(!is.na(landevap))
exist.dat.nocrop$climate_main <- substr(exist.dat.nocrop$layer_climate, 1, 1)
exist.dat.nocrop$climate_sub <- substr(exist.dat.nocrop$layer_climate, 2, 2)
exist.dat.nocrop <- exist.dat.nocrop[!duplicated(exist.dat.nocrop)]

#Remove unusable values

exist.dat<- exist.dat[is.na(HDI) == FALSE & is.na(landevap) == FALSE]
exist.dat<- exist.dat[, -c("V1", "summass", "sumarea", "mindose", "name_chemclass", "name_country")]
exist.dat$logwd <- log10(exist.dat$wd)
exist.dat <- exist.dat[, c("name_group_yy", "name_crop_earthstat", "casrn_chemical" ,"name_region.y","layer_climate","score", "HDI", "landevap", "logwd", "croparea", "sumcroparea" )]
exist.dat$climate_main <- substr(exist.dat$layer_climate, 1, 1)
exist.dat$climate_sub <- substr(exist.dat$layer_climate, 2, 2)

#filter chemicals with more scenarios, sigf = 0.05
exist.dat <- exist.dat %>%
  group_by(casrn_chemical) %>%
  mutate(count = n()) %>%
  filter(count > 60) %>%
  ungroup() %>%
  select(-count)

chem.list<- unique(exist.dat$casrn_chemical)

#Prep set of general scenarios that match chemical-specific scenarios
exist.dat.cpp <- exist.dat

general.scenario<- exist.dat.cpp %>%
  select(-casrn_chemical, -logwd) %>%
  distinct()

#adding the no crop country combination to general scenatio
general.scenario <- rbind(general.scenario, exist.dat.nocrop)
#general.scenario<- general.scenario[!duplicated(general.scenario)] #Lots of replicates in scenarios, remove

#Prep dataset to fill in with ML values
#dat.out.knnlog<- matrix(data="X", nrow=length(chem.list), ncol= 5, dimnames=list(NULL, c("CASRN", "nNo", "nYes", "rsq1", "rsq2")))
dat.out.cpp0.05 <- matrix(data = "X", nrow = length(chem.list), ncol = 9, dimnames = list(NULL, c("CASRN", "nNo", "nYes","No_ActualNo", "Yes_ActualNo", "No_ActualYes",
                                                                                               "Yes_ActualYes", "No_Accuracy", "Yes_Accuracy")))
dat.out.cpp0.1 <- matrix(data = "X", nrow = length(chem.list), ncol = 10, dimnames = list(NULL, c("CASRN", "nNo", "nYes","No_ActualNo", "Yes_ActualNo", "No_ActualYes",
                                                                                                   "Yes_ActualYes", "No_Accuracy", "Yes_Accuracy", "CPEfficiency")))
dat.out.cpp0.15 <- matrix(data = "X", nrow = length(chem.list), ncol = 10, dimnames = list(NULL, c("CASRN", "nNo", "nYes","No_ActualNo", "Yes_ActualNo", "No_ActualYes",
                                                                                                   "Yes_ActualYes", "No_Accuracy", "Yes_Accuracy", "CPEfficiency")))
dat.out.cpp0.2 <- matrix(data = "X", nrow = length(chem.list), ncol = 10, dimnames = list(NULL, c("CASRN", "nNo", "nYes","No_ActualNo", "Yes_ActualNo", "No_ActualYes",
                                                                                                   "Yes_ActualYes", "No_Accuracy", "Yes_Accuracy", "CPEfficiency")))
#dat.out.xgblog<- matrix(data="X", nrow=length(chem.list), ncol= 5, dimnames=list(NULL, c("CASRN", "nNo", "nYes", "rsq1", "rsq2")))
i<-2
pred_xgb_full <- data.frame()
no.dat.cpp.pred.exist.full0.05 <- data.frame()
no.dat.cpp.pred.exist.full0.1 <- data.frame()
no.dat.cpp.pred.exist.full0.15 <- data.frame()
no.dat.cpp.pred.exist.full0.2 <- data.frame()
chem.list <- c("133-06-2", "60-51-5", "133-06-2", "1918-02-1", "143390-89-0")
#Building a separate model for each chem - loop through chem list and test accuracy
for(i in c(1:length(chem.list))){
  chem<- chem.list[i]
  #chem <- "1314-56-3"
  exist.dat.chem <- exist.dat %>%
    filter(casrn_chemical == chem) #%>%
  present.dat <- exist.dat.chem
  present.dat<- present.dat %>%
    select(-casrn_chemical, -logwd) %>%
    distinct()
  #present.dat<- present.dat[!duplicated(present.dat)]
  absent.dat <- general.scenario %>%
    anti_join(present.dat, by = c("score", "name_crop_earthstat", "name_region.y", 
                                  "name_group_yy", "layer_climate", "landevap", "HDI", "croparea", "sumcroparea"))
  
  #Does a scenario use a chem?
  present.dat$ChemPresent<- "Yes"
  absent.dat$ChemPresent <- "No"
  
  machine.dat<- rbind(present.dat, absent.dat)
  write.csv(machine.dat, paste0("machine.dat", chem, ".csv"))
  
  #conformal prediction steps
  machine.dat.cpp <- read.csv(paste0("machine.dat", chem, ".csv"))
  machine.dat.cpp <- machine.dat.cpp %>%
    select(ChemPresent, name_crop_earthstat, name_region.y, name_group_yy, HDI, landevap, layer_climate, climate_main, climate_sub, score, croparea, sumcroparea)
    #select(ChemPresent, name_crop_earthstat, name_region.y, name_group_yy, HDI, landevap, layer_climate, climate_main, climate_sub, score, croparea, sumcroparea)
  machine.dat.cpp[,1] = as.factor(machine.dat.cpp[,1])
  machine.dat.cpp[,1] = as.numeric(unlist(machine.dat.cpp[,1]))
  
  write.csv(machine.dat.cpp, paste0("20221104machine.dat.cpp", chem, ".csv"), row.names=FALSE)
 

  no.dat$ChemPresent = "null"
  no.dat.cpp <- no.dat %>%
    select(ChemPresent, name_crop_earthstat, name_region.y, name_group_yy, HDI, landevap, layer_climate, climate_main, climate_sub, score, croparea, sumcroparea)
    #select(ChemPresent, name_crop_earthstat, name_region.y, name_group_yy, HDI, landevap, layer_climate, climate_main, climate_sub, score, croparea, sumcroparea)
  #write.csv(no.dat.cpp, "20221104no.dat.cpp.csv", row.names = FALSE)
  #Create a training set - this is what will be used to build the model (random 70% of machine data rows)
  #train.dat.cpp<- createDataPartition(machine.dat.cpp$layer_climate, p=0.7, list=FALSE)
  
  #The remaining 30% of rows are the test set - these will determine how good the model is
  #test.dat.cpp<- machine.dat.cpp[-train.dat.cpp,]
  #train.dat.cpp<- machine.dat.cpp[train.dat.cpp,]
  #scramble train dataset
  #train.dat.cpp.shuffle <- transform(train.dat.cpp, name_crop_earthstat = sample(name_crop_earthstat), HDI = sample(HDI),
                                     #landevap = sample(landevap))
  #train.dat.cpp <- train.dat.cpp.shuffle
  
  
  #pValues = ICPClassification(train.dat.cpp, test.dat.cpp,method = "rf", nrTrees = 100)
  #pValues = ICPClassification(machine.dat.cpp, no.dat.cpp, ratioTrain = 0.7, method = "rf", nrTrees = 100)
  #pValues.df <- as.data.frame(pValues)
  #testLabels = test.dat.cpp[,1]
  #CPEfficiency <- CPEfficiency(pValues, testLabels, 0.05)
  #CPErrorRate <- CPErrorRate(pValues, testLabels)
  #CPValidity <- CPValidity(pValues, testLabels)
  #CPCalibrationPlot(pValues, test.dat.cpp, "blue")
  
  #10 times conformal prediction on train/test data
  bootstrap_preds_yes <- matrix(data = "X", nrow = nrow(no.dat.cpp), ncol = 10)
  bootstrap_preds_no <- matrix(data = "X", nrow = nrow(no.dat.cpp), ncol = 10)
  for(b in c(1:10)){
    #train_idxs = sample(seq(n),size = n, replace = TRUE)
    originalData_r <- machine.dat.cpp
    #originalData_r <- originalData_r %>%
    #select(-orig.id)
    pValues = ICPClassification(originalData_r, no.dat.cpp, ratioTrain = 0.7, method = "rf", nrTrees = 100)
    for(j in c(1:nrow(no.dat.cpp))){
      bootstrap_preds_no[j,b] = pValues[j,1]
      bootstrap_preds_yes[j,b] = pValues[j,2]
      #print(i)
    }
    print(b)
  }
  
  sig <- 0.05
  bootstrap_preds_yes.df <- as.data.frame(bootstrap_preds_yes)
  bootstrap_preds_yes.df <- bootstrap_preds_yes.df %>%
    mutate_if(is.character, as.numeric) %>%
    rowwise() %>%
    mutate(medianyes = median(c_across(where(is.numeric)), na.rm = TRUE)) #%>%
    # mutate(chemPred = case_when(median >= sig ~ "V2",
    #                             median < sig ~ "V1"))
  bootstrap_preds_no.df <- as.data.frame(bootstrap_preds_no)
  bootstrap_preds_no.df <- bootstrap_preds_no.df %>%
    mutate_if(is.character, as.numeric) %>%
    rowwise() %>%
    mutate(medianno = median(c_across(where(is.numeric)), na.rm = TRUE)) #%>%
    # mutate(chemPred = case_when(median >= sig ~ "V2",
    #                             median < sig ~ "V1"))
  
  pValues.df <- cbind(bootstrap_preds_yes.df, bootstrap_preds_no.df)
  pValues.df <- pValues.df %>%
    select(medianno, medianyes) %>%
    mutate(V1 = medianno, V2 = medianyes) %>%
    select(V1, V2)
  

  

  sig <- 0.05
  pValues.df.pred0.05 <- pValues.df %>%
    mutate(chemPred = case_when(#V1 >= sig & V2 >= sig ~ "both p >= sig",
                                V1 < sig & V2 < sig ~ "both p < sig",
                                V1 >= sig & V2 < V1 ~ "V1",
                                V2 >= sig & V2 > V1 ~ "V2"))
  #chemPred <- pValues.df.pred0.05$chemPred
  #chemPred <- bootstrap_preds_yes.df$chemPred
  #test.dat.cpp$ChemPresent <- as.character(test.dat.cpp$ChemPresent)
  #cm <- as.data.table(table(test.dat.cpp$ChemPresent, chemPred))
  
  #test.dat.pred <- cbind(test.dat.cpp, pValues.df.pred0.05)
  #write.csv(test.dat.pred, "20221101test.dat.pred.csv")
  
  no.dat.cpp.pred <- cbind(no.dat.cpp, pValues.df.pred0.05)
  
  #no.dat.cpp.pred <- cbind(no.dat.cpp, bootstrap_preds_yes.df)
  #no.dat.cpp.pred$casrn_chemical <- chem
  #no.dat.cpp.pred.exist <- no.dat.cpp.pred %>%
    #filter(chemPred == "V2") %>%
    #filter(V2 >= sig) %>%
    #select(-ChemPresent, -chemPred, -V1, -V2)
  write.csv(no.dat.cpp.pred, paste0("20221104no.dat.cpp.pred0.05", chem, ".csv"), row.names = FALSE)
  #write.csv(no.dat.cpp.pred.exist, paste0("20221003no.dat.cpp.pred.exist", chem, ".csv"))
  
  #no.dat.cpp.pred.exist.full0.05 <- rbind(no.dat.cpp.pred.exist.full0.05, no.dat.cpp.pred.exist)
  # val.out <- c(chem,
  # nrow(absent.dat),
  # nrow(present.dat),
  # cm[V2 == "1" & chemPred == "V1"]$N,
  # cm[V2 == "1" & chemPred == "V2"]$N,
  # cm[V2 == "2" & chemPred == "V1"]$N,
  # cm[V2 == "2" & chemPred == "V2"]$N,
  # cm[V2 == "1" & chemPred == "V1"]$N/length(which(test.dat.cpp$ChemPresent == "1")),
  # cm[V2 == "2" & chemPred == "V2"]$N/length(which(test.dat.cpp$ChemPresent == "2")))
  # 
  # write.csv(cm, "20221101cpp_cm.csv")

  #dat.out.cpp0.05[i,] <- val.out
  
  # sig <- 0.1
  # pValues.df.pred0.1 <- pValues.df %>%
  #   mutate(chemPred = case_when(V1 >= sig & V2 >= sig ~ "both p >= sig",
  #                               V1 < sig & V2 < sig ~ "both p < sig",
  #                               V1 >= sig & V2 < V1 ~ "V1",
  #                               V2 >= sig & V2 > V1 ~ "V2"))
  #chemPred <- pValues.df.pred0.05$chemPred
  #chemPred <- bootstrap_preds_yes.df$chemPred
  #test.dat.cpp$ChemPresent <- as.character(test.dat.cpp$ChemPresent)
  #cm <- as.data.table(table(test.dat.cpp$ChemPresent, chemPred))
  
  
  # no.dat.cpp.pred <- cbind(no.dat.cpp, pValues.df.pred0.1)
  # #no.dat.cpp.pred <- cbind(no.dat.cpp, bootstrap_preds_yes.df)
  # no.dat.cpp.pred$casrn_chemical <- chem
  # no.dat.cpp.pred.exist <- no.dat.cpp.pred %>%
  #   filter(chemPred == "V2") %>%
  #   #filter(V2 >= sig) %>%
  #   select(-ChemPresent, -chemPred, -V1, -V2)
  #write.csv(no.dat.cpp.pred, paste0("20221004no.dat.cpp.pred0.05", chem, ".csv"))
  #write.csv(no.dat.cpp.pred.exist, paste0("20221003no.dat.cpp.pred.exist", chem, ".csv"))
  
  #no.dat.cpp.pred.exist.full0.1 <- rbind(no.dat.cpp.pred.exist.full0.1, no.dat.cpp.pred.exist)
  # # val.out <- c(chem,
  # #              nrow(absent.dat),
  # #              nrow(present.dat),
  # #              cm[V2 == "1" & chemPred == "V1"]$N,
  # #              cm[V2 == "1" & chemPred == "V2"]$N,
  # #              cm[V2 == "2" & chemPred == "V1"]$N,
  # #              cm[V2 == "2" & chemPred == "V2"]$N,
  # #              cm[V2 == "1" & chemPred == "V1"]$N/length(which(test.dat.cpp$ChemPresent == "1")),
  # #              cm[V2 == "2" & chemPred == "V2"]$N/length(which(test.dat.cpp$ChemPresent == "2")),
  # #              CPEfficiency)
  # # 
  # # dat.out.cpp0.1[i,] <- val.out
  # 
  # sig <- 0.15
  # pValues.df.pred0.15 <- pValues.df %>%
  #   mutate(chemPred = case_when(V1 >= sig & V2 >= sig ~ "both p >= sig",
  #                               V1 < sig & V2 < sig ~ "both p < sig",
  #                               V1 >= sig & V2 < V1 ~ "V1",
  #                               V2 >= sig & V2 > V1 ~ "V2"))
  # chemPred <- pValues.df.pred0.15$chemPred
  # #test.dat.cpp$ChemPresent <- as.character(test.dat.cpp$ChemPresent)
  # #cm <- as.data.table(table(test.dat.cpp$ChemPresent, chemPred))
  # 
  # no.dat.cpp.pred <- cbind(no.dat.cpp, pValues.df.pred0.15)
  # no.dat.cpp.pred$casrn_chemical <- chem
  # no.dat.cpp.pred.exist <- no.dat.cpp.pred %>%
  #   #filter(chemPred == "V2") %>%
  #   filter(V2 >= sig) %>%
  #   select(-ChemPresent, -chemPred, -V1, -V2)
  # #write.csv(no.dat.cpp.pred, paste0("20221004no.dat.cpp.pred0.15", chem, ".csv"))
  # #write.csv(no.dat.cpp.pred.exist, paste0("20221003no.dat.cpp.pred.exist", chem, ".csv"))
  # 
  # no.dat.cpp.pred.exist.full0.15 <- rbind(no.dat.cpp.pred.exist.full0.15, no.dat.cpp.pred.exist)
  # # val.out <- c(chem,
  # #              nrow(absent.dat),
  # #              nrow(present.dat),
  # #              cm[V2 == "1" & chemPred == "V1"]$N,
  # #              cm[V2 == "1" & chemPred == "V2"]$N,
  # #              cm[V2 == "2" & chemPred == "V1"]$N,
  # #              cm[V2 == "2" & chemPred == "V2"]$N,
  # #              cm[V2 == "1" & chemPred == "V1"]$N/length(which(test.dat.cpp$ChemPresent == "1")),
  # #              cm[V2 == "2" & chemPred == "V2"]$N/length(which(test.dat.cpp$ChemPresent == "2")),
  # #              CPEfficiency)
  # # 
  # # dat.out.cpp0.15[i,] <- val.out
  # 
  # sig <- 0.2
  # pValues.df.pred0.2 <- pValues.df %>%
  #   mutate(chemPred = case_when(V1 >= sig & V2 >= sig ~ "both p >= sig",
  #                               V1 < sig & V2 < sig ~ "both p < sig",
  #                               V1 >= sig & V2 < V1 ~ "V1",
  #                               V2 >= sig & V2 > V1 ~ "V2"))
  # chemPred <- pValues.df.pred0.2$chemPred
  # #test.dat.cpp$ChemPresent <- as.character(test.dat.cpp$ChemPresent)
  # #cm <- as.data.table(table(test.dat.cpp$ChemPresent, chemPred))
  # 
  # no.dat.cpp.pred <- cbind(no.dat.cpp, pValues.df.pred0.2)
  # no.dat.cpp.pred$casrn_chemical <- chem
  # no.dat.cpp.pred.exist <- no.dat.cpp.pred %>%
  #   #filter(chemPred == "V2") %>%
  #   filter(V2 >= sig) %>%
  #   select(-ChemPresent, -chemPred, -V1, -V2)
  # #write.csv(no.dat.cpp.pred, paste0("20221004no.dat.cpp.pred0.2", chem, ".csv"))
  # #write.csv(no.dat.cpp.pred.exist, paste0("20221003no.dat.cpp.pred.exist", chem, ".csv"))
  # 
  # no.dat.cpp.pred.exist.full0.2 <- rbind(no.dat.cpp.pred.exist.full0.2, no.dat.cpp.pred.exist)
  # val.out <- c(chem,
  #              nrow(absent.dat),
  #              nrow(present.dat),
  #              cm[V2 == "1" & chemPred == "V1"]$N,
  #              cm[V2 == "1" & chemPred == "V2"]$N,
  #              cm[V2 == "2" & chemPred == "V1"]$N,
  #              cm[V2 == "2" & chemPred == "V2"]$N,
  #              cm[V2 == "1" & chemPred == "V1"]$N/length(which(test.dat.cpp$ChemPresent == "1")),
  #              cm[V2 == "2" & chemPred == "V2"]$N/length(which(test.dat.cpp$ChemPresent == "2")),
  #              CPEfficiency)

  #dat.out.cpp0.2[i,] <- val.out
  
  exist.dat.xgb <- exist.dat.chem 
  machine.dat.xgb <- exist.dat.xgb %>%
     select(-casrn_chemical)
  #no.dat.xgb <- no.dat.cpp.pred.exist #%>%
     #select(-casrn_chemical)
  
  write.csv(machine.dat.xgb, paste0("20221104machine.dat.xgb.", chem,"csv"), row.names = FALSE)
  
  #Create a training set - this is what will be used to build the model (random 70% of machine data rows)
  #train.dat.xgb1<- createDataPartition(machine.dat.xgb$layer_climate, p=0.7, list=FALSE)
  
  #The remaining 30% of rows are the test set - these will determine how good the model is
  #test.dat.xgb<- machine.dat.xgb[-train.dat.xgb1,] 
  #train.dat.xgb<- machine.dat.xgb[train.dat.xgb1,]
   
  # control<- trainControl(method="repeatedcv", number=5, repeats=2) #This will be the control for each model run
  # #fit.model<- train(ChemPresent ~., data=train.dat, method="xgbTree", metric="Accuracy", trControl=control)
  # #fit.model <- train(logwd ~., data = train.dat.xgb, method = "xgbTree", metric = "RMSE", trControl = control)
  # fit.model <- train(logwd ~., data = machine.dat.xgb, method = "xgbTree", metric = "RMSE", trControl = control)
  # # #xgb.results <- fit.model$results
  # # 
  # #predict for test data
  # #newdata <- test.dat.xgb %>%
  #   #select(-logwd)
  # #newdata$name_crop_earthstat[which(!(newdata$name_crop_earthstat %in% unique(train.dat.xgb$name_crop_earthstat)))] <- NA
  # #newdata$name_region[which(!(newdata$name_region.y %in% unique(train.dat.xgb$name_region.y)))] <- NA
  # 
  # newdata.xgb <- no.dat.xgb
  # newdata.xgb$name_crop_earthstat[which(!(newdata.xgb$name_crop_earthstat %in% unique(machine.dat.xgb$name_crop_earthstat)))] <- NA
  # newdata.xgb$name_region.y[which(!(newdata.xgb$name_region.y %in% unique(machine.dat.xgb$name_region.y)))] <- NA
  # newdata.xgb$layer_climate[which(!(newdata.xgb$layer_climate %in% unique(machine.dat.xgb$layer_climate)))] <- NA
  # newdata.xgb$climate_sub[which(!(newdata.xgb$climate_sub %in% unique(machine.dat.xgb$climate_sub)))] <- NA
  # newdata.xgb$climate_main[which(!(newdata.xgb$climate_main %in% unique(machine.dat.xgb$climate_main)))] <- NA
  # newdata.xgb$name_group_yy[which(!(newdata.xgb$name_group_yy %in% unique(machine.dat.xgb$name_group_yy)))] <- NA
  # newdata <- newdata.xgb %>%
  #   filter(!is.na(.))
  # pred.xgb <- predict(fit.model, newdata = newdata, na.action = na.pass)
  # 
  # chem_pred_xgb_full <- cbind(newdata, pred.xgb)
  # 
  # chem_pred_xgb_full$casrn_chemical <- chem
  # 
  # chem_pred_xgb_full <- cbind(newdata, test.dat.xgb$logwd, chem_pred_xgb_full$pred)
  # 
  # write.csv(chem_pred_xgb_full, paste0("20221101chem_pred_xgb_full_no", chem, ".csv"))
  # 
  # pred_xgb_full <- rbind(pred_xgb_full, chem_pred_xgb_full)
  # 
  # 
  # pred_orig <- as.data.frame(cbind(test.dat.xgb$logwd, chem_pred_xgb_full$pred))
  # pred_orig %>%
  #   ggplot(aes(x = 10^(test.dat.xgb$logwd), y = 10^(chem_pred_xgb_full$pred))) +
  #   geom_point() +
  #   geom_abline() +
  #   scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #                 labels = trans_format("log10", math_format(10^.x)))+
  #   scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #                 labels = trans_format("log10", math_format(10^.x)))+
  #   labs(title = paste0("pred_orig", chem))
  # ggsave(paste0("20220914pred_orig", chem, ".jpeg"))
  
  print(i)
}

fit.model

write.csv(no.dat.cpp.pred.exist.full, "20220915no.dat.cpp.pred.exist.full.csv")
write.csv(dat.out.cpp0.05, "20221028dat.out.cpp0.05_withoutnorows.csv")
write.csv(no.dat.cpp.pred.exist.full0.1, "20221022no.dat.cpp.pred.exist.nocfropareasig0.1_2.csv")

dat.out.cpp0.05_withnorows <- read.csv("20221027dat.out.cpp0.05_withnorows.csv")
dat.out.cpp0.05_withnorows <- as.data.frame(dat.out.cpp0.05_withnorows)
dat.out.cpp0.05_withnorows$source <- "predict_withextrarows"
dat.out.cpp0.05_withnorows <- dat.out.cpp0.05_withnorows %>%
  select(-X)
dat.out.cpp0.05_withoutnorows <- read.csv("20221028dat.out.cpp0.05_withoutnorows.csv")
dat.out.cpp0.05_withoutnorows <- dat.out.cpp0.05_withoutnorows %>%
  select(-X)
dat.out.cpp0.05_scramble <- read.csv("20221027dat.out.cpp0.05_withnorowsscramble.csv")
dat.out.cpp0.05_scramble$source <- "scramble"
dat.out.cpp0.05_scramble <- dat.out.cpp0.05_scramble %>%
  select(-X)
dat.out.cpp0.05_withoutnorows$source <- "predict_withoutextrarows"
dat.out.cpp0.05_all <- rbind(dat.out.cpp0.05_withnorows, dat.out.cpp0.05_scramble, dat.out.cpp0.05_withoutnorows)
dat.out.cpp0.05_all %>%
  ggplot(aes(x = as.numeric(Yes_Accuracy), y = as.numeric(No_Accuracy))) +
  geom_point(aes(color = source), alpha = 0.5)

dat.out.cppover0.8 <- dat.out.cpp0.05_all %>%
  filter(Yes_Accuracy > 0.8 & No_Accuracy > 0.8)
ggsave("20221028scramblepredictaccuracy.jpeg", height = 20, width = 20, units = "cm")

dat.out.cpp <- dat.out.cpp0.05 %>%
  left_join(dat.out.cpp0.1, by = c("CASRN")) %>%
  left_join(dat.out.cpp0.15, by = c("CASRN")) %>%
  left_join(dat.out.cpp0.2, by = c("CASRN"))
              


dat.out.cpp0.2 <- as.data.frame(dat.out.cpp0.2)
dat.out.cpp %>%
  ggplot() +
  geom_point(aes(x = as.numeric(Yes_Accuracy.x), y = as.numeric(No_Accuracy.x)), color = "blue") +
  geom_point(aes(x = as.numeric(Yes_Accuracy.y), y = as.numeric(No_Accuracy.y)), color = "red") +
  geom_point(aes(x = as.numeric(Yes_Accuracy.x.x), y = as.numeric(No_Accuracy.x.x)), color = "yellow") +
  geom_point(aes(x = as.numeric(Yes_Accuracy.y.y), y = as.numeric(No_Accuracy.y.y)), color = "orange") 
write.csv(dat.out.cpp, "20220912dat.out.cpp_ccclevel.csv")
write.csv(dat.out.xgblog, "20220816dat.out.xgblog_cropspecific.csv")


#predict dose for all chems (no separation for treatment methods)
exist.dat<- fread("C:/me/data gaps/20220912country_input_full_all.csv")
#Remove unusable values
exist.dat<- exist.dat[is.na(HDI) == FALSE & is.na(landevap) == FALSE]
exist.dat<- exist.dat[, -c("V1", "summass", "sumarea", "mindose", "name_chemclass", "name_country")]
exist.dat$logwd <- log10(exist.dat$wd)
exist.dat <- exist.dat[, c("name_cropgroup", "name_crop_earthstat", "casrn_chemical" ,"name_region.y","layer_climate","score", "HDI", "landevap", "logwd" )]
exist.dat$climate_main <- substr(exist.dat$layer_climate, 1, 1)
exist.dat$climate_sub <- substr(exist.dat$layer_climate, 2, 2)


#filter chemicals with more scenarios, sigf = 0.05
exist.dat <- exist.dat %>%
  group_by(casrn_chemical) %>%
  mutate(count = n()) %>%
  filter(count > 10) %>%
  ungroup() %>%
  select(-count)

chem.list<- unique(exist.dat$casrn_chemical)

dat.out.xgblog<- matrix(data="X", nrow=length(chem.list), ncol= 5, dimnames=list(NULL, c("CASRN", "nNo", "nYes", "rsq1", "rsq2")))

i<-303
for(i in c(648:length(chem.list))){
  chem<- chem.list[i]
  exist.dat.xgb <- exist.dat %>%
    filter(casrn_chemical == chem) #%>%
  
  #only do regression model on existing chemical scenarios
  machine.dat.xgb <- exist.dat.xgb %>%
    select(-casrn_chemical)

  
  #Create a training set - this is what will be used to build the model (random 70% of machine data rows)
  train.dat.xgb1<- createDataPartition(machine.dat.xgb$layer_climate, p=0.7, list=FALSE)
  
  #The remaining 30% of rows are the test set - these will determine how good the model is
  test.dat.xgb<- machine.dat.xgb[-train.dat.xgb1,] 
  train.dat.xgb<- machine.dat.xgb[train.dat.xgb1,]
  
  #xgboost to predict using caret
  control<- trainControl(method="repeatedcv", number=5, repeats=2) #This will be the control for each model run
  #fit.model<- train(ChemPresent ~., data=train.dat, method="xgbTree", metric="Accuracy", trControl=control)
  #fit.model <- train(logwd ~., data = train.dat.xgb, method = "xgbTree", metric = "RMSE", trControl = control)
  fit.model <- train(logwd ~., data = train.dat.xgb, method = "xgbTree", metric = "RMSE", trControl = control)
  #xgb.results <- fit.model$results
  
  #predict for test data
  newdata <- test.dat.xgb %>%
    select(-logwd)
  newdata$name_crop_earthstat[which(!(newdata$name_crop_earthstat %in% unique(train.dat.xgb$name_crop_earthstat)))] <- NA
  newdata$name_region.y[which(!(newdata$name_region.y %in% unique(train.dat.xgb$name_region.y)))] <- NA
  newdata$layer_climate[which(!(newdata$layer_climate %in% unique(train.dat.xgb$layer_climate)))] <- NA
  newdata <- newdata %>%
    filter(!is.na(.))
  
  pred.xgb <- predict(fit.model, newdata = newdata, na.action = na.pass)
  
  chem_pred_xgb_full <- cbind(newdata, pred.xgb)
  
  chem_pred_xgb_full$casrn_chemical <- chem
  
  write.csv(chem_pred_xgb_full, paste0("20220912chem_pred_xgb_full", chem, ".csv"))
  
  pred <- chem_pred_xgb_full$pred
  orig <- test.dat.xgb$logwd

  rss <- sum((pred - orig)^2)
  tss <- sum((orig-mean(orig))^2)
  rsq1 <- 1 - rss/tss
  rsq2 <- cor(pred, orig)
  
  pred_orig <- as.data.frame(cbind(test.dat.xgb$logwd, chem_pred_xgb_full$pred))
  pred_orig %>%
    ggplot(aes(x = 10^(test.dat.xgb$logwd), y = 10^(chem_pred_xgb_full$pred))) +
    geom_point() +
    geom_abline() +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))+
    labs(title = paste0("pred_orig", chem))
  ggsave(paste0("20220914pred_orig", chem, ".jpeg"))
  
  val.out<- c(chem,
              nrow(exist.dat), #Total Nos in dataset
              nrow(exist.dat.xgb), #Total yes in dataset
              rsq1,
              rsq2)

  dat.out.xgblog[i,]<- val.out
  
  print(i)
}

write.csv(dat.out.xgblog,"20220915pred_orig_xgboost.csv")

# mass dose prediction
#Existed: complete info including chem use for crop-country combos
exist.dat<- fread("C:/me/data gaps/20220926country_input_full_all_massdose.csv")
exist.dat <- exist.dat %>%
  left_join(crop_id, by = c("name_crop_earthstat" = "CROPNAME"))
no.country.continent <- fread("C:/me/Database/harmonization/no.country.continent.csv")
no.dat <- fread("C:/me/data gaps/20220926country_input_n_all_massdose.csv")
no.dat <- no.dat %>%
  filter(score != 0) %>%
  left_join(crop_id, by = c("variable" = "CROPNAME")) %>%
  mutate(name_group_yy = name_group_yy)%>%
  left_join(no.country.continent, by = c("name_country" = "name_country")) %>%
  mutate(name_region.y = name_region) %>%
  mutate(name_crop_earthstat = variable) %>%
  select("name_group_yy", "name_crop_earthstat", "name_region.y","layer_climate","score", "HDI", "landevap" ) %>%
  filter(!is.na(HDI)) %>%
  filter(!is.na(landevap))
no.dat$climate_main <- substr(no.dat$layer_climate, 1, 1)
no.dat$climate_sub <- substr(no.dat$layer_climate, 2, 2)
no.dat <- no.dat[!duplicated(no.dat)]


#Remove unusable values

exist.dat<- exist.dat[is.na(HDI) == FALSE & is.na(landevap) == FALSE]
exist.dat$logwd <- log10(exist.dat$wd)
exist.dat <- exist.dat[, c("name_group_yy", "name_crop_earthstat", "casrn_chemical" ,"name_region.y","layer_climate","score", "HDI", "landevap", "logwd" )]
exist.dat$climate_main <- substr(exist.dat$layer_climate, 1, 1)
exist.dat$climate_sub <- substr(exist.dat$layer_climate, 2, 2)

#filter chemicals with more scenarios, sigf = 0.05
exist.dat <- exist.dat %>%
  group_by(casrn_chemical) %>%
  mutate(count = n()) %>%
  filter(count > 60) %>%
  ungroup() %>%
  select(-count)

chem.list<- unique(exist.dat$casrn_chemical)

#Prep set of general scenarios that match chemical-specific scenarios
exist.dat.cpp <- exist.dat

general.scenario<- exist.dat.cpp %>%
  select(-casrn_chemical, -logwd) %>%
  distinct()
#general.scenario<- general.scenario[!duplicated(general.scenario)] #Lots of replicates in scenarios, remove

#Prep dataset to fill in with ML values
#dat.out.knnlog<- matrix(data="X", nrow=length(chem.list), ncol= 5, dimnames=list(NULL, c("CASRN", "nNo", "nYes", "rsq1", "rsq2")))
dat.out.cpp <- matrix(data = "X", nrow = length(chem.list), ncol = 11, dimnames = list(NULL, c("CASRN", "nNo", "nYes","No_ActualNo", "Yes_ActualNo", "No_ActualYes",
                                                                                               "Yes_ActualYes", "No_Accuracy", "Yes_Accuracy", "CPEfficiency", "CPErrorRate")))
dat.out.xgblog<- matrix(data="X", nrow=length(chem.list), ncol= 5, dimnames=list(NULL, c("CASRN", "nNo", "nYes", "rsq1", "rsq2")))
i<-3
pred_xgb_full <- data.frame()
no.dat.cpp.pred.exist.full <- data.frame()
#Building a separate model for each chem - loop through chem list and test accuracy
for(i in c(641:length(chem.list))){
  chem <- "1071-83-6"
  #chem<- chem.list[i]
  exist.dat.chem <- exist.dat %>%
    filter(casrn_chemical == chem)
  present.dat <- exist.dat.chem
  present.dat<- present.dat %>%
    select(-casrn_chemical, -logwd) %>%
    distinct()
  #present.dat<- present.dat[!duplicated(present.dat)]
  absent.dat <- general.scenario %>%
    anti_join(present.dat, by = c("score", "name_crop_earthstat", "name_region.y", 
                                  "name_group_yy", "layer_climate", "landevap", "HDI"))
  
  #Does a scenario use a chem?
  present.dat$ChemPresent<- "Yes"
  absent.dat$ChemPresent <- "No"
  
  machine.dat<- rbind(present.dat, absent.dat)
  write.csv(machine.dat, paste0("machine.dat", chem, ".csv"))
  
  #conformal prediction steps
  machine.dat.cpp <- read.csv(paste0("machine.dat", chem, ".csv"))
  machine.dat.cpp <- machine.dat.cpp %>%
    select(ChemPresent, name_crop_earthstat, name_region.y, name_group_yy, HDI, landevap, layer_climate, climate_main, climate_sub, score)
  machine.dat.cpp[,1] = as.factor(machine.dat.cpp[,1])
  machine.dat.cpp[,1] = as.numeric(unlist(machine.dat.cpp[,1]))
  
  no.dat$ChemPresent = "null"
  no.dat.cpp <- no.dat %>%
    select(ChemPresent, name_crop_earthstat, name_region.y, name_group_yy, HDI, landevap, layer_climate, climate_main, climate_sub, score)
  #Create a training set - this is what will be used to build the model (random 70% of machine data rows)
  #train.dat.cpp1<- createDataPartition(machine.dat.cpp$name_region.y, p=0.7, list=FALSE)
  
  #The remaining 30% of rows are the test set - these will determine how good the model is
  #test.dat.cpp<- machine.dat.cpp[-train.dat.cpp1,]
  #train.dat.cpp<- machine.dat.cpp[train.dat.cpp1,]
  
  #pValues = ICPClassification(train.dat.cpp, test.dat.cpp,method = "rf", nrTrees = 100)
  pValues = ICPClassification(machine.dat.cpp, no.dat.cpp,method = "rf", nrTrees = 100)
  pValues.df <- as.data.frame(pValues)
  #testLabels = test.dat.cpp[,1]
  #CPEfficiency <- CPEfficiency(pValues, testLabels, 0.05)
  #CPErrorRate <- CPErrorRate(pValues, testLabels)
  #CPValidity <- CPValidity(pValues, testLabels)
  
  sig <- 0.05
  pValues.df.pred <- pValues.df %>%
    mutate(chemPred = case_when(V1 >= sig & V2 >= sig ~ "both p >= sig",
                                V1 < sig & V2 < sig ~ "both p < sig",
                                V1 >= sig & V2 < V1 ~ "V1",
                                V2 >= sig & V2 > V1 ~ "V2"))
  #chemPred <- pValues.df.pred$chemPred
  #test.dat.cpp$ChemPresent <- as.character(test.dat.cpp$ChemPresent)
  #cm <- as.data.table(table(test.dat.cpp$ChemPresent, chemPred))
  #no.data.cpp.pred <- cbind(test.dat.cpp, pValues.df.pred)
  #no.data.cpp.pred$casrn_chemical <- chem
  # no.data.cpp.pred.exist <- no.data.cpp.pred %>%
  #   filter(chemPred == "V2") %>%
  #   select(-ChemPresent, -chemPred, -V1, -V2)
  no.dat.cpp.pred <- cbind(no.dat.cpp, pValues.df.pred)
  no.dat.cpp.pred$casrn_chemical <- chem
  no.dat.cpp.pred.exist <- no.dat.cpp.pred %>%
    filter(chemPred == "V2") %>%
    select(-ChemPresent, -chemPred, -V1, -V2)
  #write.csv(no.data.cpp.pred, paste0("20220921no.dat.cpp.pred", chem, ".csv"))
  #write.csv(no.data.cpp.pred.exist, paste0("20220921no.data.cpp.pred.exist", chem, ".csv"))
  
  #no.dat.cpp.pred.exist.full <- rbind(no.dat.cpp.pred.exist.full, no.dat.cpp.pred.exist)
  # val.out <- c(chem,
  # nrow(absent.dat),
  # nrow(present.dat),
  # cm[V2 == "1" & chemPred == "V1"]$N,
  # cm[V2 == "1" & chemPred == "V2"]$N,
  # cm[V2 == "2" & chemPred == "V1"]$N,
  # cm[V2 == "2" & chemPred == "V2"]$N,
  # cm[V2 == "1" & chemPred == "V1"]$N/length(which(test.dat.cpp$ChemPresent == "1")),
  # cm[V2 == "2" & chemPred == "V2"]$N/length(which(test.dat.cpp$ChemPresent == "2")),
  # CPEfficiency, CPErrorRate)
  # 
  # dat.out.cpp[i,] <- val.out
  
  machine.dat.xgb <- exist.dat.chem %>%
    select(-casrn_chemical)
  #Create a training set - this is what will be used to build the model (random 70% of machine data rows)
  #train.dat.xgb1<- createDataPartition(machine.dat.xgb$name_region.y, p=0.7, list=FALSE)
  
  #The remaining 30% of rows are the test set - these will determine how good the model is
  #test.dat.xgb<- machine.dat.xgb[-train.dat.xgb1,] 
  #train.dat.xgb<- machine.dat.xgb[train.dat.xgb1,]
  #machine.dat.xgb <- exist.dat.xgb %>%
    #filter(ChemPresent == "Yes") %>%
    #select(-ChemPresent, -summass)
  no.dat.xgb <- no.dat.cpp.pred.exist %>%
    select(-casrn_chemical) 
  #no.dat.xgb <- test.dat.xgb %>%
    #filter(ChemPresent == "Yes") %>%
    #select(-ChemPresent, -logwd)
  #machine.dat.xgb <- train.dat.xgb
  control<- trainControl(method="repeatedcv", number=5, repeats=2) #This will be the control for each model run
  #fit.model<- train(ChemPresent ~., data=train.dat, method="xgbTree", metric="Accuracy", trControl=control)
  #fit.model <- train(logwd ~., data = train.dat.xgb, method = "xgbTree", metric = "RMSE", trControl = control)
  fit.model <- train(logwd ~., data = machine.dat.xgb, method = "xgbTree", metric = "RMSE", trControl = control)
  #xgb.results <- fit.model$results
  
  #predict for test data
  #newdata <- test.dat.xgb %>%
  #select(-logwd)
  #newdata$name_crop_earthstat[which(!(newdata$name_crop_earthstat %in% unique(train.dat.xgb$name_crop_earthstat)))] <- NA
  #newdata$name_region[which(!(newdata$name_region.y %in% unique(train.dat.xgb$name_region.y)))] <- NA
  
  #newdata.xgb <- no.dat.xgb
  newdata.xgb <- no.dat.xgb
  newdata.xgb$name_crop_earthstat[which(!(newdata.xgb$name_crop_earthstat %in% unique(machine.dat.xgb$name_crop_earthstat)))] <- NA
  newdata.xgb$name_region.y[which(!(newdata.xgb$name_region.y %in% unique(machine.dat.xgb$name_region.y)))] <- NA
  newdata.xgb$name_group_yy[which(!(newdata.xgb$name_group_yy %in% unique(machine.dat.xgb$name_group_yy)))] <- NA
  newdata.xgb$layer_climate[which(!(newdata.xgb$layer_climate %in% unique(machine.dat.xgb$layer_climate)))] <- NA
  newdata.xgb$climate_sub[which(!(newdata.xgb$climate_sub %in% unique(machine.dat.xgb$climate_sub)))] <- NA
  newdata.xgb$climate_main[which(!(newdata.xgb$climate_main %in% unique(machine.dat.xgb$climate_main)))] <- NA
  
  newdata <- newdata.xgb %>%
    filter(!is.na(.))
  #newdata_input <- newdata %>%
    #select(-logwd)
  pred.xgb <- predict(fit.model, newdata = newdata, na.action = na.pass)
  
  chem_pred_xgb_full <- cbind(newdata, pred.xgb)
  
  
  # pred <- chem_pred_xgb_full$pred.xgb
  # orig <- chem_pred_xgb_full$logwd
  # 
  # rss <- sum((pred - orig)^2)
  # tss <- sum((orig-mean(orig))^2)
  # rsq1 <- 1 - rss/tss
  # rsq2 <- cor(pred, orig)
  # 
  # 
  # chem_pred_xgb_full %>%
  #   ggplot(aes(x = 10^(pred.xgb), y = 10^(logwd))) +
  #   geom_point() +
  #   geom_abline() +
  #   scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #                 labels = trans_format("log10", math_format(10^.x)))+
  #   scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #                 labels = trans_format("log10", math_format(10^.x)))+
  #   labs(title = paste0("pred_orig", chem))
  # ggsave(paste0("20220926pred_orig", chem, ".jpeg"))
  
  # val.out<- c(chem,
  #             nrow(exist.dat.chem), #Total Nos in dataset
  #             nrow(test.dat.xgb), #Total yes in dataset
  #             rsq1,
  #             rsq2)
  
  #dat.out.xgblog[i,]<- val.out
  
  chem_pred_xgb_full$casrn_chemical <- chem
  
  write.csv(chem_pred_xgb_full, paste0("20221008chem_pred_xgb_full", chem, ".csv"))
  
  pred_xgb_full <- rbind(pred_xgb_full, chem_pred_xgb_full)
  
  
  print(i)
}



# USA mass dose prediction
#Existed: complete info including chem use for crop-country combos
exist.dat<- fread("C:/me/data gaps/20221028country_input_full_usa.csv")
#exist.dat <- fread("C:/me/data gaps/20220929country_input_full_all_massdose_usa.csv")
exist.dat.nocrop <- fread("C:/me/data gaps/20221028country_imputation_combi_existcountrynocrop_morocco.csv")
exist.dat <- exist.dat %>%
  select(-name_group_yy) %>%
  mutate(name_region.y = case_when(name_country == "CENTRAL AMERICA-CARIBBEAN" ~ "North America",
                                   TRUE ~ name_region.y)) %>%
  left_join(crop_id, by = c("name_crop_earthstat" = "CROPNAME")) %>%
  left_join(country_croparea, by = c("name_country", "name_crop_earthstat" = "variable")) %>%
  left_join(no.country.continent, by = c("name_country"))
no.country.continent <- read_xlsx("C:/me/Database/harmonization/no.country.continent.xlsx", sheet = 1)
no.dat <- fread("C:/me/data gaps/20221028country_input_n_usa.csv")
#no.dat <- fread("C:/me/data gaps/20220929country_input_n_all_massdose_usa.csv")
no.dat <- no.dat %>%
  filter(score != 0) %>%
  mutate(name_country = name_cntr) %>%
  #filter(name_country %in% country_noprediction$name_country) %>%
  left_join(crop_id, by = c("variable" = "CROPNAME")) %>%
  mutate(name_group_yy = name_group_yy)%>%
  left_join(no.country.continent, by = c("name_country" = "name_country")) %>%
  left_join(country_croparea, by = c("name_country", "variable")) %>%
  mutate(name_region.y = name_region) %>%
  mutate(name_crop_earthstat = variable) %>%
  select("name_group_yy", "name_crop_earthstat", "name_region.y","layer_climate","score", "HDI", "landevap", "croparea", "sumcroparea", "FDI_in", "GCF" ) %>%
  filter(!is.na(HDI)) %>%
  filter(!is.na(landevap))
no.dat$climate_main <- substr(no.dat$layer_climate, 1, 1)
no.dat$climate_sub <- substr(no.dat$layer_climate, 2, 2)
no.dat <- no.dat[!duplicated(no.dat)]
exist.dat.nocrop <- exist.dat.nocrop %>%
  filter(score != 0) %>%
  mutate(name_country = name_cntr) %>%
  select("layer_climate", "name_country", "score", "variable", "HDI", "landevap") %>%
  left_join(crop_id, by = c("variable" = "CROPNAME")) %>%
  mutate(name_group_yy = name_group_yy)%>%
  left_join(no.country.continent, by = c("name_country" = "name_country")) %>%
  left_join(country_croparea, by = c("name_country", "variable")) %>%
  mutate(name_region.y = name_region) %>%
  mutate(name_crop_earthstat = variable) %>%
  select("name_group_yy", "name_crop_earthstat", "name_region.y","layer_climate","score", "HDI", "landevap", "croparea", "sumcroparea" ) %>%
  filter(!is.na(HDI)) %>%
  filter(!is.na(landevap))
exist.dat.nocrop$climate_main <- substr(exist.dat.nocrop$layer_climate, 1, 1)
exist.dat.nocrop$climate_sub <- substr(exist.dat.nocrop$layer_climate, 2, 2)
exist.dat.nocrop <- exist.dat.nocrop[!duplicated(exist.dat.nocrop)]


#Remove unusable values

exist.dat<- exist.dat[is.na(HDI) == FALSE & is.na(landevap) == FALSE]
exist.dat$logwd <- log10(exist.dat$wd)
exist.dat <- exist.dat[, c("name_group_yy", "name_crop_earthstat", "casrn_chemical" ,"name_region.y","layer_climate","score", "HDI", "landevap", "logwd", "croparea", "sumcroparea", "FDI_in", "GCF")]
exist.dat$climate_main <- substr(exist.dat$layer_climate, 1, 1)
exist.dat$climate_sub <- substr(exist.dat$layer_climate, 2, 2)

#filter chemicals with more scenarios, sigf = 0.05
exist.dat <- exist.dat %>%
  group_by(casrn_chemical) %>%
  mutate(count = n()) %>%
  filter(count > 60) %>%
  ungroup() %>%
  select(-count)

chem.list<- unique(exist.dat$casrn_chemical)

#Prep set of general scenarios that match chemical-specific scenarios
exist.dat.cpp <- exist.dat

general.scenario<- exist.dat.cpp %>%
  select(-casrn_chemical, -logwd) %>%
  distinct()
#general.scenario<- general.scenario[!duplicated(general.scenario)] #Lots of replicates in scenarios, remove

#adding the no crop country combination to general scenatio
general.scenario <- rbind(general.scenario, exist.dat.nocrop)

#Prep dataset to fill in with ML values
#dat.out.knnlog<- matrix(data="X", nrow=length(chem.list), ncol= 5, dimnames=list(NULL, c("CASRN", "nNo", "nYes", "rsq1", "rsq2")))
dat.out.cpp <- matrix(data = "X", nrow = length(chem.list), ncol = 11, dimnames = list(NULL, c("CASRN", "nNo", "nYes","No_ActualNo", "Yes_ActualNo", "No_ActualYes",
                                                                                               "Yes_ActualYes", "No_Accuracy", "Yes_Accuracy", "CPEfficiency", "CPErrorRate")))
dat.out.xgblog<- matrix(data="X", nrow=length(chem.list), ncol= 5, dimnames=list(NULL, c("CASRN", "nNo", "nYes", "rsq1", "rsq2")))
i<-3
pred_xgb_full <- data.frame()
no.dat.cpp.pred.exist.full0.05 <- data.frame()
no.dat.cpp.pred.exist.full0.1 <- data.frame()
no.dat.cpp.pred.exist.full0.15 <- data.frame()
no.dat.cpp.pred.exist.full0.2 <- data.frame()
#Building a separate model for each chem - loop through chem list and test accuracy
for(i in c(1:length(chem.list))){
  chem<- chem.list[i]
  #chem <- "1314-56-3"
  exist.dat.chem <- exist.dat %>%
    filter(casrn_chemical == chem)
  present.dat <- exist.dat.chem
  present.dat<- present.dat %>%
    select(-casrn_chemical, -logwd) %>%
    distinct()
  #present.dat<- present.dat[!duplicated(present.dat)]
  absent.dat <- general.scenario %>%
    anti_join(present.dat, by = c("score", "name_crop_earthstat", "name_region.y", 
                                  "name_group_yy", "layer_climate", "landevap", "HDI", "croparea", "sumcroparea", "FDI_in", "GCF"))
  
  #Does a scenario use a chem?
  present.dat$ChemPresent<- "Yes"
  absent.dat$ChemPresent <- "No"
  
  machine.dat<- rbind(present.dat, absent.dat)
  write.csv(machine.dat, paste0("machine.dat", chem, ".csv"))
  
  #conformal prediction steps
  machine.dat.cpp <- read.csv(paste0("machine.dat", chem, ".csv"))
  machine.dat.cpp <- machine.dat.cpp %>%
    select(ChemPresent, name_crop_earthstat, name_region.y, name_group_yy, HDI, landevap, layer_climate, climate_main, climate_sub, score, croparea, sumcroparea, FDI_in, GCF)
  machine.dat.cpp[,1] = as.factor(machine.dat.cpp[,1])
  machine.dat.cpp[,1] = as.numeric(unlist(machine.dat.cpp[,1]))
  
  no.dat$ChemPresent = "null"
  no.dat.cpp <- no.dat %>%
    select(ChemPresent, name_crop_earthstat, name_region.y, name_group_yy, HDI, landevap, layer_climate, climate_main, climate_sub, score, croparea, sumcroparea, FDI_in, GCF)
  #Create a training set - this is what will be used to build the model (random 70% of machine data rows)
  #train.dat.cpp1<- createDataPartition(machine.dat.cpp$name_region.y, p=0.7, list=FALSE)
  
  write.csv(machine.dat.cpp, paste0("machine.dat.cpp_E",chem, ".csv"))
  write.csv(no.dat.cpp, paste0("no.dat.cpp_E", chem, ".csv"))
  #The remaining 30% of rows are the test set - these will determine how good the model is
  #test.dat.cpp<- machine.dat.cpp[-train.dat.cpp1,]
  #train.dat.cpp<- machine.dat.cpp[train.dat.cpp1,]
  
  #pValues = ICPClassification(train.dat.cpp, test.dat.cpp,method = "rf", nrTrees = 100)
  pValues = ICPClassification(machine.dat.cpp, no.dat.cpp,method = "rf", nrTrees = 100)
  pValues.df <- as.data.frame(pValues)
  #testLabels = test.dat.cpp[,1]
  #CPEfficiency <- CPEfficiency(pValues, testLabels, 0.05)
  #CPErrorRate <- CPErrorRate(pValues, testLabels)
  #CPValidity <- CPValidity(pValues, testLabels)
  
  sig <- 0.05
  pValues.df.pred <- pValues.df %>%
    mutate(chemPred = case_when(#V1 >= sig & V2 >= sig ~ "both p >= sig",
                                V1 < sig & V2 < sig ~ "both p < sig",
                                V1 >= sig & V2 < V1 ~ "V1",
                                V2 >= sig & V2 > V1 ~ "V2"))
  #chemPred <- pValues.df.pred$chemPred
  #test.dat.cpp$ChemPresent <- as.character(test.dat.cpp$ChemPresent)
  #cm <- as.data.table(table(test.dat.cpp$ChemPresent, chemPred))
  #no.data.cpp.pred <- cbind(test.dat.cpp, pValues.df.pred)
  #no.data.cpp.pred$casrn_chemical <- chem
  # no.data.cpp.pred.exist <- no.data.cpp.pred %>%
  #   filter(chemPred == "V2") %>%
  #   select(-ChemPresent, -chemPred, -V1, -V2)
  no.dat.cpp.pred <- cbind(no.dat.cpp, pValues.df.pred)
  no.dat.cpp.pred$casrn_chemical <- chem
  no.dat.cpp.pred.exist <- no.dat.cpp.pred %>%
    filter(chemPred == "V2") %>%
    select(-ChemPresent, -chemPred)
  #write.csv(no.data.cpp.pred, paste0("20220921no.dat.cpp.pred", chem, ".csv"))
  #write.csv(no.dat.cpp.pred.exist, paste0("20220929no.data.cpp.pred.exist", chem, ".csv"))
  no.dat.cpp.pred.exist.full0.05 <- rbind(no.dat.cpp.pred.exist.full0.05, no.dat.cpp.pred.exist)
  
  
  # sig <- 0.1
  # pValues.df.pred <- pValues.df %>%
  #   mutate(chemPred = case_when(V1 >= sig & V2 >= sig ~ "both p >= sig",
  #                               V1 < sig & V2 < sig ~ "both p < sig",
  #                               V1 >= sig & V2 < V1 ~ "V1",
  #                               V2 >= sig & V2 > V1 ~ "V2"))
  # #chemPred <- pValues.df.pred$chemPred
  # #test.dat.cpp$ChemPresent <- as.character(test.dat.cpp$ChemPresent)
  # #cm <- as.data.table(table(test.dat.cpp$ChemPresent, chemPred))
  # #no.data.cpp.pred <- cbind(test.dat.cpp, pValues.df.pred)
  # #no.data.cpp.pred$casrn_chemical <- chem
  # # no.data.cpp.pred.exist <- no.data.cpp.pred %>%
  # #   filter(chemPred == "V2") %>%
  # #   select(-ChemPresent, -chemPred, -V1, -V2)
  # no.dat.cpp.pred <- cbind(no.dat.cpp, pValues.df.pred)
  # no.dat.cpp.pred$casrn_chemical <- chem
  # no.dat.cpp.pred.exist <- no.dat.cpp.pred %>%
  #   filter(chemPred == "V2") %>%
  #   select(-ChemPresent, -chemPred, -V1, -V2)
  # #write.csv(no.data.cpp.pred, paste0("20220921no.dat.cpp.pred", chem, ".csv"))
  # #write.csv(no.dat.cpp.pred.exist, paste0("20220929no.data.cpp.pred.exist", chem, ".csv"))
  # no.dat.cpp.pred.exist.full0.1 <- rbind(no.dat.cpp.pred.exist.full0.1, no.dat.cpp.pred.exist)
  # 
  # sig <- 0.15
  # pValues.df.pred <- pValues.df %>%
  #   mutate(chemPred = case_when(#V1 >= sig & V2 >= sig ~ "both p >= sig",
  #                               V1 < sig & V2 < sig ~ "both p < sig",
  #                               V1 >= sig & V2 < V1 ~ "V1",
  #                               V2 >= sig & V2 > V1 ~ "V2"))
  # #chemPred <- pValues.df.pred$chemPred
  # #test.dat.cpp$ChemPresent <- as.character(test.dat.cpp$ChemPresent)
  # #cm <- as.data.table(table(test.dat.cpp$ChemPresent, chemPred))
  # #no.data.cpp.pred <- cbind(test.dat.cpp, pValues.df.pred)
  # #no.data.cpp.pred$casrn_chemical <- chem
  # # no.data.cpp.pred.exist <- no.data.cpp.pred %>%
  # #   filter(chemPred == "V2") %>%
  # #   select(-ChemPresent, -chemPred, -V1, -V2)
  # no.dat.cpp.pred <- cbind(no.dat.cpp, pValues.df.pred)
  # no.dat.cpp.pred$casrn_chemical <- chem
  # no.dat.cpp.pred.exist <- no.dat.cpp.pred %>%
  #   filter(chemPred == "V2") %>%
  #   select(-ChemPresent, -chemPred, -V1, -V2)
  # #write.csv(no.data.cpp.pred, paste0("20220921no.dat.cpp.pred", chem, ".csv"))
  # #write.csv(no.dat.cpp.pred.exist, paste0("20220929no.data.cpp.pred.exist", chem, ".csv"))
  # no.dat.cpp.pred.exist.full0.15 <- rbind(no.dat.cpp.pred.exist.full0.15, no.dat.cpp.pred.exist)
  # 
  # sig <- 0.2
  # pValues.df.pred <- pValues.df %>%
  #   mutate(chemPred = case_when(V1 >= sig & V2 >= sig ~ "both p >= sig",
  #                               V1 < sig & V2 < sig ~ "both p < sig",
  #                               V1 >= sig & V2 < V1 ~ "V1",
  #                               V2 >= sig & V2 > V1 ~ "V2"))
  # #chemPred <- pValues.df.pred$chemPred
  # #test.dat.cpp$ChemPresent <- as.character(test.dat.cpp$ChemPresent)
  # #cm <- as.data.table(table(test.dat.cpp$ChemPresent, chemPred))
  # #no.data.cpp.pred <- cbind(test.dat.cpp, pValues.df.pred)
  # #no.data.cpp.pred$casrn_chemical <- chem
  # # no.data.cpp.pred.exist <- no.data.cpp.pred %>%
  # #   filter(chemPred == "V2") %>%
  # #   select(-ChemPresent, -chemPred, -V1, -V2)
  # no.dat.cpp.pred <- cbind(no.dat.cpp, pValues.df.pred)
  # no.dat.cpp.pred$casrn_chemical <- chem
  # no.dat.cpp.pred.exist <- no.dat.cpp.pred %>%
  #   filter(chemPred == "V2") %>%
  #   select(-ChemPresent, -chemPred, -V1, -V2)
  # #write.csv(no.data.cpp.pred, paste0("20220921no.dat.cpp.pred", chem, ".csv"))
  # #write.csv(no.dat.cpp.pred.exist, paste0("20220929no.data.cpp.pred.exist", chem, ".csv"))
  # no.dat.cpp.pred.exist.full0.2 <- rbind(no.dat.cpp.pred.exist.full0.2, no.dat.cpp.pred.exist)
  # machine.dat.xgb <- exist.dat.chem %>%
  #   select(-casrn_chemical)
  # #Create a training set - this is what will be used to build the model (random 70% of machine data rows)
  # #train.dat.xgb1<- createDataPartition(machine.dat.xgb$name_region.y, p=0.7, list=FALSE)
  # 
  # #The remaining 30% of rows are the test set - these will determine how good the model is
  # #test.dat.xgb<- machine.dat.xgb[-train.dat.xgb1,] 
  # #train.dat.xgb<- machine.dat.xgb[train.dat.xgb1,]
  # #machine.dat.xgb <- exist.dat.xgb %>%
  # #filter(ChemPresent == "Yes") %>%
  # #select(-ChemPresent, -summass)
  # no.dat.xgb <- no.dat.cpp.pred.exist %>%
  #   select(-casrn_chemical) 
  # #no.dat.xgb <- test.dat.xgb %>%
  # #filter(ChemPresent == "Yes") %>%
  # #select(-ChemPresent, -logwd)
  # #machine.dat.xgb <- train.dat.xgb
  # control<- trainControl(method="repeatedcv", number=5, repeats=2) #This will be the control for each model run
  # #fit.model<- train(ChemPresent ~., data=train.dat, method="xgbTree", metric="Accuracy", trControl=control)
  # #fit.model <- train(logwd ~., data = train.dat.xgb, method = "xgbTree", metric = "RMSE", trControl = control)
  # fit.model <- train(logwd ~., data = machine.dat.xgb, method = "xgbTree", metric = "RMSE", trControl = control)
  # #xgb.results <- fit.model$results
  # 
  # #predict for test data
  # #newdata <- test.dat.xgb %>%
  # #select(-logwd)
  # #newdata$name_crop_earthstat[which(!(newdata$name_crop_earthstat %in% unique(train.dat.xgb$name_crop_earthstat)))] <- NA
  # #newdata$name_region[which(!(newdata$name_region.y %in% unique(train.dat.xgb$name_region.y)))] <- NA
  # 
  # #newdata.xgb <- no.dat.xgb
  # newdata.xgb <- no.dat.xgb
  # newdata.xgb$name_crop_earthstat[which(!(newdata.xgb$name_crop_earthstat %in% unique(machine.dat.xgb$name_crop_earthstat)))] <- NA
  # newdata.xgb$name_region.y[which(!(newdata.xgb$name_region.y %in% unique(machine.dat.xgb$name_region.y)))] <- NA
  # newdata.xgb$name_group_yy[which(!(newdata.xgb$name_group_yy %in% unique(machine.dat.xgb$name_group_yy)))] <- NA
  # newdata.xgb$layer_climate[which(!(newdata.xgb$layer_climate %in% unique(machine.dat.xgb$layer_climate)))] <- NA
  # newdata.xgb$climate_sub[which(!(newdata.xgb$climate_sub %in% unique(machine.dat.xgb$climate_sub)))] <- NA
  # newdata.xgb$climate_main[which(!(newdata.xgb$climate_main %in% unique(machine.dat.xgb$climate_main)))] <- NA
  # 
  # newdata <- newdata.xgb %>%
  #   filter(!is.na(.))
  # #newdata_input <- newdata %>%
  # #select(-logwd)
  # pred.xgb <- predict(fit.model, newdata = newdata, na.action = na.pass)
  # 
  # chem_pred_xgb_full <- cbind(newdata, pred.xgb)
  # 
  # 
  # # pred <- chem_pred_xgb_full$pred.xgb
  # # orig <- chem_pred_xgb_full$logwd
  # # 
  # # rss <- sum((pred - orig)^2)
  # # tss <- sum((orig-mean(orig))^2)
  # # rsq1 <- 1 - rss/tss
  # # rsq2 <- cor(pred, orig)
  # # 
  # # 
  # # chem_pred_xgb_full %>%
  # #   ggplot(aes(x = 10^(pred.xgb), y = 10^(logwd))) +
  # #   geom_point() +
  # #   geom_abline() +
  # #   scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
  # #                 labels = trans_format("log10", math_format(10^.x)))+
  # #   scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  # #                 labels = trans_format("log10", math_format(10^.x)))+
  # #   labs(title = paste0("pred_orig", chem))
  # # ggsave(paste0("20220926pred_orig", chem, ".jpeg"))
  # 
  # # val.out<- c(chem,
  # #             nrow(exist.dat.chem), #Total Nos in dataset
  # #             nrow(test.dat.xgb), #Total yes in dataset
  # #             rsq1,
  # #             rsq2)
  # 
  # #dat.out.xgblog[i,]<- val.out
  # 
  # chem_pred_xgb_full$casrn_chemical <- chem
  # 
  # write.csv(chem_pred_xgb_full, paste0("20220929chem_pred_xgb_full_usa", chem, ".csv"))
  # 
  # pred_xgb_full <- rbind(pred_xgb_full, chem_pred_xgb_full)
  
  print(i)
}




unique(no.dat.cpp.pred.exist.full0.05$casrn_chemical)

write.csv(no.dat.cpp.pred.exist.full0.05, "20221103no.dat.cpp.pred.exist.full0.05_m.csv")
usa_crop_chemcount <- no.dat.cpp.pred.exist.full0.2 %>%
  group_by(name_crop_earthstat) %>%
  summarize(countchem = n_distinct(casrn_chemical))
usa_crop_chemcount %>%
  ggplot(aes(x = countchem)) +
  geom_histogram()
usa_crop_chemcount <- input_caschem %>%
  filter(name_country == "United States of America") %>%
  group_by(name_crop_earthstat) %>%
  summarize(countchem = n_distinct(casrn_chemical))
