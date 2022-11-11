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

setwd("C:/Users/markos/Documents/DataDownloads/combi_cntr_crop_pest")
#Existed: complete info including chem use for crop-country combos in EU 
exist.dat<- fread("C:/me/data gaps/20220725country_input_full.csv")
no.dat <- fread("C:/me/data gaps/20220725country_input_n.csv")

head(exist.dat)
summary(exist.dat)

table(exist.dat$name_crop)
table(exist.dat$variable)
table(exist.dat$name_cntr)
table(exist.dat$name_ai)
table(exist.dat$layer_climate)
table(exist.dat$regu_eu)
table(exist.dat$regu_cntr)

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

exist.dat<- exist.dat[is.na(HDI) == FALSE & is.na(landevap) == FALSE]
exist.dat<- exist.dat[, -c("V1", "summass", "sumarea", "mindose", "name_chemclass", "name_country", "name_crop_earthstat")]
exist.dat$climate_main <- substr(exist.dat$layer_climate, 1, 1)
exist.dat$climate_sub <- substr(exist.dat$layer_climate, 2, 2)

exist.dat.chem <- exist.dat %>%
  group_by(casrn_chemical) %>%
  mutate(count = n()) %>%
  filter(count > 60)

#Remove unusable values
length(unique(exist.dat$name_ai)) == length(unique(exist.dat$casrn_chemical)) #Only use 1 of these columns

chem.list<- unique(exist.dat.chem$casrn_chemical)

#Prep set of general scenarios that match chemical-specific scenarios
exist.dat.cpp <- exist.dat[,-c("wd")]
general.scenario<- exist.dat.cpp[, -c("casrn_chemical")]
general.scenario<- general.scenario[!duplicated(general.scenario)] #Lots of replicates in scenarios, remove

#Prep dataset to fill in with ML values
dat.out.knnlog<- matrix(data="X", nrow=length(chem.list), ncol= 5, dimnames=list(NULL, c("CASRN", "nNo", "nYes", "rsq1", "rsq2")))
dat.out.cpp <- matrix(data = "X", nrow = length(chem.list), ncol = 10, dimnames = list(NULL, c("CASRN", "nNo", "nYes","No_ActualNo", "Yes_ActualNo", "No_ActualYes",
                                                                                               "Yes_ActualYes", "No_Accuracy", "Yes_Accuracy", "CPEfficiency")))
dat.out<- matrix(data="X", nrow=length(chem.list), ncol= 5, dimnames=list(NULL, c("CASRN", "nNo", "nYes", "rsq1", "rsq2")))
i<-1
#Building a separate model for each chem - loop through chem list and test accuracy
#for(i in c(1:length(chem.list))){
for(i in c(1:length(chem.list))){
  chem<- chem.list[i]
  exist.dat.chem <- exist.dat %>%
    filter(casrn_chemical == chem) #%>%
    #filter(casrn_chemical == "51-03-6")
  
  #present.dat<- exist.dat.cpp[casrn_chemical %in% chem] #ID scenarios that use that chem
  present.dat <- exist.dat.chem
  present.dat<- present.dat[, -c("casrn_chemical")]
  present.dat<- present.dat[!duplicated(present.dat)]
  #x<- row.match(present.dat, general.scenario) #ID scenarios that do not use that chem
  #absent.dat<- general.scenario[!na.omit(x)]
  absent.dat <- general.scenario %>%
    anti_join(present.dat, by = c("score", "Short Name", "name_region.y", 
                                  "name_cropgroup", "layer_climate", "landevap", "HDI"))
  #Do any of the present rows match absent rows? If yes, there's an error
  #test<- row.match(present.dat, absent.dat)
  #if(length(which(is.na(test) == FALSE)) != 0) stop ('Multi-value match')
  
  #Does a scenario use a chem?
  present.dat$ChemPresent<- "Yes"
  absent.dat$ChemPresent <- "No"
  absent.dat$wd<- 0
  
  #machine.dat<- as.data.table(rbind(present.dat, absent.dat))
  machine.dat <- rbind(present.dat, absent.dat)
  
  #Set machine dat types to match necessary input for ML 
  # machine.dat<- transform(
  #   machine.dat,
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
  #   #name_country=as.factor(name_country),
  #   score = as.numeric(score),
  #   name_crop_earthstat=as.factor(name_crop_earthstat), 
  #   name_region = as.factor(name_region),
  #   name_cropgroup = as.factor(name_cropgroup),
  #   #wd = as.numeric(wd),
  #   ChemPresent = as.factor(ChemPresent)
  # )
  
  write.csv(machine.dat, paste0("machine.dat", chem, ".csv"))
  machine.dat.cpp <- read.csv("machine.dat.csv")
  machine.dat.cpp <- machine.dat.cpp %>%
    select(ChemPresent, name_crop_earthstat, name_region, name_cropgroup, HDI, landevap, layer_climate, climate_main, climate_sub, score)
  machine.dat.cpp[,1] = as.factor(machine.dat.cpp[,1])
  machine.dat.cpp[,1] = as.numeric(unlist(machine.dat.cpp[,1]))
  
  #Create a training set - this is what will be used to build the model (random 70% of machine data rows)
  train.dat.cpp<- createDataPartition(machine.dat.cpp$layer_climate, p=0.7, list=FALSE)
  
  #The remaining 30% of rows are the test set - these will determine how good the model is
  test.dat.cpp<- machine.dat.cpp[-train.dat.cpp,] 
  
  train.dat.cpp<- machine.dat.cpp[train.dat.cpp,] 
  
  pValues = ICPClassification(train.dat.cpp, test.dat.cpp,method = "rf", nrTrees = 100)
  pValues.df <- as.data.frame(pValues)
  testLabels = test.dat[,1]
  CPEfficiency <- CPEfficiency(pValues, testLabels, 0.05)
  #CPErrorRate <- CPErrorRate(pValues, testLabels)
  #CPValidity <- CPValidity(pValues, testLabels)
  
  sig <- 0.05
  pValues.df.pred <- pValues.df %>%
    mutate(chemPred = case_when(V1 >= sig & V2 >= sig ~ "both p >= sig",
                                V1 < sig & V2 < sig ~ "both p < sig",
                                V1 >= sig & V2 < V1 ~ "V1",
                                V2 >= sig & V2 > V1 ~ "V2"))
  chemPred <- pValues.df.pred$chemPred
  test.dat.cpp$ChemPresent <- as.character(test.dat.cpp$ChemPresent)
  cm <- as.data.table(table(test.dat.cpp$ChemPresent, chemPred))
  
  val.out <- c(chem,
               nrow(absent.dat),
               nrow(present.dat),
               cm[V2 == "1" & chemPred == "V1"]$N,
               cm[V2 == "1" & chemPred == "V2"]$N,
               cm[V2 == "2" & chemPred == "V1"]$N,
               cm[V2 == "2" & chemPred == "V2"]$N,
               cm[V2 == "1" & chemPred == "V1"]$N/length(which(test.dat.cpp$ChemPresent == "1")),
               cm[V2 == "2" & chemPred == "V2"]$N/length(which(test.dat.cpp$ChemPresent == "2")),
               CPEfficiency)
  
  dat.out.cpp[i,] <- val.out
  
  #test.dat.cpp.pred <- cbind(test.dat.cpp, chemPred)
  
  #machine.dat.xgb <- machine.dat %>%
    #filter(ChemPresent == "Yes")
  
  machine.dat.xgb <- present.dat
  machine.dat.xgb <- machine.dat.xgb %>%
    select(-ChemPresent) %>%
    mutate(logwd = log10(wd)) %>%
    select(-wd)
  
  #Create a training set - this is what will be used to build the model (random 70% of machine data rows)
  train.dat.xgb1<- createDataPartition(machine.dat.xgb$layer_climate, p=0.7, list=FALSE)
  
  #The remaining 30% of rows are the test set - these will determine how good the model is
  test.dat.xgb<- machine.dat.xgb[-train.dat.xgb1,] 
  train.dat.xgb<- machine.dat.xgb[train.dat.xgb1,]
  

  dat.xgb.rec <- 
    recipe(logwd ~., data = data_pred_train) %>%
    step_dummy(all_nominal())
  
  dat.xgb.boot <- data_pred_train %>%
    bootstraps()
  
  #dat.xgb.boot
  
  dat.xgb <- 
    boost_tree() %>%
    set_mode("regression") %>%
    set_engine("xgboost")
  
  #define a basic xgboost model
  fit_bootstrap <- function(index){
    xgb_boot <- dat.xgb.boot$splits[[index]] %>%
      training()
    
    workflow() %>%
      add_recipe(dat.xgb.rec) %>%
      add_model(dat.xgb) %>%
      #boost_tree() %>%
      #set_mode("regression") %>%
      #set_engine("xgboost") %>%
      fit(xgb_boot) %>%
      write_rds(paste0("model_", index, ".rds"))
  }
  
  #fit to 25 bootstrapped datasets
  for(i in 1:25) {
    fit_bootstrap(i)
  }
  
  predict_bootstrap <- function(new_data, index) {
    read_rds(paste0("model_", index, ".rds")) %>%
      predict(new_data) %>%
      rename(!!sym(paste0("pred_", index)) := .pred)
  }
  
  training_preds <- 
    seq(1,25) %>%
    map_dfc(~predict_bootstrap(data_pred_train, .x))
  
  #training_preds
  
  training_preds %>%
    bind_cols(data_pred_train) %>%
    rowid_to_column() %>%
    pivot_longer(starts_with("pred_"),
                 names_to = "model",
                 values_to = ".pred") %>%
    group_by(rowid) %>%
    summarise(logwd = max(logwd),
              pred_mean = mean(.pred),
              std_dev = sd(.pred)) %>%
    riekelib::normal_interval(pred_mean, std_dev) %>%
    ggplot(aes(x = logwd,
               y = pred_mean)) +
    geom_point(alpha = 0.5) +
    geom_segment(aes(x = logwd,
                 xend = logwd,
                 y = ci_lower,
                 yend = ci_upper),
                 alpha = 0.25)
  
  
 p<-  seq(1,25) %>%
    map_dfc(~predict_bootstrap(data_pred_test, .x)) %>%
    bind_cols(data_pred_test) %>%
    rowid_to_column() %>%
    pivot_longer(starts_with("pred_"),
                 names_to = "model",
                 values_to = ".pred") %>%
    group_by(rowid) %>%
    summarise(logwd = max(logwd),
              pred_mean = mean(.pred),
              std_dev = sd(.pred)) %>%
    riekelib::normal_interval(pred_mean, std_dev) %>%
    ggplot(aes(x = logwd,
               y= pred_mean)) +
    geom_point(alpha = 0.5) +
    geom_segment(aes(x = logwd,
                     xend = logwd,
                     y = ci_lower,
                     yend = ci_upper),
                 alpha = 0.25)
  p <- seq(1,25) %>%
    map_dfc(~predict_bootstrap(data_pred_test, .x)) %>%
    bind_cols(data_pred_test) %>%
    rowid_to_column() %>%
    pivot_longer(starts_with("pred_"),
                 names_to = "model",
                 values_to = ".pred") %>%
    group_by(rowid) %>%
    summarise(logwd = max(logwd),
              pred_mean = mean(.pred),
              std_dev = sd(.pred)) %>%
    riekelib::normal_interval(pred_mean, std_dev)
    
  
  pred_2 <- p$pred_mean
  orig_2 <- p$logwd
  rss <- sum((pred_2 - orig_2)^2)
  tss <- sum((orig_2-mean(orig_2))^2)
  rsq1 <- 1 - rss/tss
  rsq2 <- cor(p$logwd, p$pred_mean)
  
  control<- trainControl(method="repeatedcv", number=5, repeats=2) #This will be the control for each model run
  #fit.model<- train(ChemPresent ~., data=train.dat, method="xgbTree", metric="Accuracy", trControl=control)
  fit.model <- train(logwd ~., data = train.dat.xgb, method = "xgbTree", metric = "RMSE", trControl = control)
  xgb.results <- fit.model$results
  
  newdata <- test.dat.xgb %>%
    select(-logwd)
  newdata$name_crop_earthstat[which(!(newdata$name_crop_earthstat %in% unique(train.dat.xgb$name_crop_earthstat)))] <- NA
  newdata$name_region[which(!(newdata$name_region %in% unique(train.dat.xgb$name_region)))] <- NA
  pred <- predict(fit.model, newdata = newdata)
  #pred_interval <- predict(fit.model$finalModel, newdata = newdata, interval = "confidence")
  
  #bootstraping
 
  
  
  #pred.df <- as.data.table(test.dat$wd, pred)
  pred.df <- cbind(test.dat.xgb, pred)
  pred.df <- as.data.frame(pred.df)
  pred.df %>%
    ggplot(aes(x = logwd, y = pred)) +
    geom_point() + 
    geom_abline() +
    ggtitle(paste0("xgb_log", chem))
  ggsave(paste0("20220630pred_dflog", chem,".png"))
  
  pred_2 <- pred.df$pred
  orig_2 <- pred.df$logwd
  rss <- sum((pred_2 - orig_2)^2)
  tss <- sum((orig_2-mean(orig_2))^2)
  rsq1 <- 1 - rss/tss
  rsq2 <- cor(pred.df$logwd, pred.df$pred)
  
  write.csv(pred.df, paste0("20220630pred_dffull_yeschemlog",chem,".csv"))
  #How well was the test set predicted? develop confusion matrix
  #cm <- as.data.table(table(test.dat$ChemPresent, pred))
  #cm <- as.data.table(table(test.dat$wd, pred))
  
  #Put together data values to track accuracy of models
  val.out<- c(chem,
              nrow(absent.dat), #Total Nos in dataset
              nrow(present.dat), #Total yes in dataset
              rsq1,
              rsq2
  ) 
  
  dat.out[i,]<- val.out
  
  #knn
  #machine.dat<- as.data.table(rbind(present.dat, absent.dat))
  machine.dat.knn <- present.dat
  
  # put outcome in its own object
  #machine.dat$ChemPresent <- as.factor(machine.dat$ChemPresent)
  #chem.outcome <- machine.dat %>% 
  #select(ChemPresent)
  
  #only do regression model on existing chemical scenarios
  machine.dat.knn$logwd <- as.numeric(log10(machine.dat.knn$wd))
  machine.dat.knn <- machine.dat.knn %>%
    select(-wd)
  #machine.dat$wd <- as.numeric(machine.dat$wd)
  chem.outcome <- machine.dat.knn %>%
    select(logwd)
  
  # remove original variable from the data set
  machine.dat.knn <- machine.dat.knn %>% 
    select(-logwd, - ChemPresent)
  
  #str(machine.dat)
  machine.dat.knn[, c("landevap")] <- scale(machine.dat.knn[, c("landevap")])
  machine.dat.knn[,c("HDI")] <- scale(machine.dat.knn[, c("HDI")])
  machine.dat.knn[,c("score")] <- scale(machine.dat.knn[, c("score")])
  #head(machine.dat)
  
  #str(machine.dat)
  
  name_crop_earthstat <- as.data.frame(dummy.code(machine.dat.knn$name_crop_earthstat))
  name_region <- as.data.frame(dummy.code(machine.dat.knn$name_region))
  name_cropgroup <- as.data.frame(dummy.code(machine.dat.knn$name_cropgroup))
  layer_climate <- as.data.frame(dummy.code(machine.dat.knn$layer_climate))
  #ChemPresent <- as.data.frame(dummy.code(machine.dat$ChemPresent))
  #climate_main <- as.data.frame(dummy.code(machine.dat$climate_main))
  #climate_sub <- as.data.frame(dummy.code(machine.dat$name_climate_sub))
  
  #without crop_earthstat
  #machine.dat <- cbind(machine.dat, name_crop_earthstat, name_region, name_cropgroup,
  #layer_climate, climate_main)
  machine.dat.knn <- cbind(machine.dat.knn, name_crop_earthstat,name_region, name_cropgroup,
                       layer_climate)
  
  machine.dat.knn <- machine.dat.knn %>% 
    select(-one_of(c("name_crop_earthstat","name_region", "name_cropgroup",
                     "layer_climate", "climate_main", "climate_sub")))
  
  class_pred_train <- machine.dat.knn[train.dat.xgb1, ]
  label_pred_train <- as.matrix(chem.outcome[train.dat.xgb1])
  data_pred_train <- as.matrix(cbind(class_pred_train, label_pred_train))
  
  class_pred_test <- machine.dat.knn[-train.dat.xgb1, ]
  label_pred_test <- as.matrix(chem.outcome[-train.dat.xgb1])
  data_pred_test <- as.matrix(cbind(class_pred_test, label_pred_test))
  
  chem_outcome_train <- chem.outcome[train.dat.xgb1, ]
  chem_outcome_test <- chem.outcome[-train.dat.xgb1, ]
  
  #chem_pred_knn <- FNN::knn.reg(train = class_pred_train, 
                                #test = class_pred_test, y = chem_outcome_train$logwd, k=20)
  dtrain.xgb <- xgb.DMatrix(data = data_pred_train, label = label_pred_train)
  dtest.xgb <- xgb.DMatrix(data = data_pred_test, label = label_pred_test)
  model <- xgboost(data = dtrain.xgb, 
                   nround = 25,
                   objective = "reg:linear")
  predictions <- predict(model, dtest.xgb)
  plot(predictions, label_pred_test)
  cor.test(predictions, label_pred_test)
  #Use the model to predict the test set\]]
  #pred <- predict(rf, newdata=test.dat[,-"ChemPresent"])
  #pred <- predict(fit.model, newdata=test.dat[,-"ChemPresent"])
  
  
  #Predict_knn <- as.data.frame(chem_pred_knn$pred)
  
  #confidence interval
  chem_pred_caret <- train(class_pred_train, chem_outcome_train$logwd, 
                           method = "knn", preProcess = c("center","scale"))
  #chem_pred_caret
  #chem_pred_caret$finalModel
  pred_interval <- predict(chem_pred_caret$finalModel, newdata = class_pred_test, interval = "confidence")
  pred <- predict(chem_pred_caret, newdata = class_pred_test)
  pred_interval_knn <- as.data.frame(pred_interval)
  pred_knn <- as.data.frame(pred)
  pred_knn_interval <- cbind(pred_knn, pred_interval_knn)
  pred_knn_interval_full <- cbind(class_pred_test, chem_outcome_test, pred_knn_interval)
  
  
  
  write.csv(pred_knn_interval, paste0("20220704pred_knn_interval", chem, ".csv"))
  
  chem_pred_knn_full<-cbind(class_pred_test, chem_outcome_test, pred_knn)
  
  #chem_pred_knn_full <- cbind(class_pred_test, chem_outcome_test,Predict_knn)
  
  chem_pred_knn_full %>%
    ggplot(aes(x = logwd, y = pred)) +
    geom_point() +
    geom_abline() +
    ggtitle(paste0("knn_log", chem))
  
  ggsave(paste0("20220704pred_df_knn_chempresentlog", chem, ".png"))
  
  write.csv(chem_pred_knn_full, paste0("20220704pred_df_knn_chempresentlog", chem, ".csv"))
  
  pred <- chem_pred_knn_full$pred
  orig <- chem_pred_knn_full$logwd
  
  rss <- sum((pred - orig)^2)
  tss <- sum((orig-mean(orig))^2)
  rsq1 <- 1 - rss/tss
  rsq2 <- cor(pred, orig)
  
  
  #knnPredict <- predict(chem_pred_knn, newdata = class_pred_test) 
  
  #confusionMatrix(knnPredict, chem_outcome_test$ChemPresent)
  
  #cm <- as.data.table(table(as.character(chem_outcome_test$ChemPresent), as.character(knnPredict)))
  
  val.out<- c(chem,
              nrow(absent.dat), #Total Nos in dataset
              nrow(present.dat), #Total yes in dataset
              rsq1,
              rsq2
  ) 
  
  dat.out.knnlog[i,]<- val.out
  print(i)
}




dat.out.xgb <- as.data.table(dat.out)
dat.out.cpp <- as.data.table(dat.out.cpp)

write.csv(dat.out.cpp, "20220630dat.out.cpp.csv")
write.csv(dat.out.xgb, "20220630dat.out.xgb.csv")

dat.out.full <- dat.out.xgb %>%
  filter(CASRN != "X") %>%
  left_join(dat.out.knnlog, by = c("CASRN"))
dat.out.knnlog <- as.data.table(dat.out.knnlog)

dat.out.full.select <- dat.out.full %>%
  mutate(rsq.xgb = round(rsq2.x,2), rsq.knn = round(rsq2.y,2)) %>%
  select(CASRN, nNo.x, nYes.x, rsq.xgb, rsq.knn) %>%
  mutate()

dat.out.full.knnbetter <- dat.out.full %>%
  filter(rsq2.y > rsq2.x)

dat.out.full.knnbetter.select <- dat.out.full.knnbetter %>%
  mutate(rsq.xgb = rsq2.x, rsq.knn = rsq2.y) %>%
  select(CASRN, nNo.x, nYes.x, rsq.xgb, rsq.knn) 

dat.out.full.knnlarger <- dat.out.full %>%
  filter(rsq.y > rsq.x)

write.csv(dat.out.full.knnbetter, "20220630dat.out.full.knnbetter.csv")

dat.out.cppoutput1 <- as.data.table(dat.out.cpp)
dat.output1 <- as.data.table(dat.out)
#write.csv(machine.dat, "chem101463-69-8.csv")
pred_2 <- pred.df$pred
orig_2 <- pred.df$wd
rss <- sum((pred_2 - orig_2)^2)
tss <- sum((orig_2-mean(orig_2))^2)
rsq <- 1 - rss/tss

pred.df3 <- as.data.table(test.dat$wd, pred)
pred.df3 <- cbind(test.dat, pred)
pred.df3 <- as.data.frame(pred.df1)
pred.df3 %>%
  ggplot(aes(x = wd, y = pred)) +
  geom_point() + 
  geom_abline()
ggsave("20220611pred_df3.png")
pred.df6 <- pred.df6 %>%
  mutate(ratio = as.numeric(pred)/as.numeric(V1))
pred.df6 %>%
  ggplot(aes(x = ratio)) +
  geom_histogram() +
  scale_x_log10()

write.csv(pred.df2, "20220610pred_df2full.csv")
write.csv(pred.df6, "pred_df6.csv")


dat.out<- as.data.table(dat.out)
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
fwrite(as.data.table(dat.out, keep.rownames = TRUE), file="211215-CM_country2.csv")
