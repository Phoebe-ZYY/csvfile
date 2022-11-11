############################################################################

### Chem imputation
###
### Author: Yuyue Zhang

############################################################################

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

#Existed: complete info including chem use for crop-country combos
exist.dat<- fread("C:/me/data gaps/20220824country_input_full.csv")
#exist.dat<- fread("C:/me/data gaps/20220908country_input_fullsoil.csv")
no.country.continent <- fread("C:/me/Database/harmonization/no.country.continent.csv")
no.dat <- fread("C:/me/data gaps/20220824country_input_n.csv")
#no.dat <- fread("C:/me/data gaps/20220908country_input_nsoil.csv")
no.dat <- no.dat %>%
  filter(score != 0) %>%
  left_join(crop_id, by = c("variable" = "CROPNAME")) %>%
  mutate(name_group_yy = name_group_yy.y)%>%
  left_join(no.country.continent, by = c("name_country" = "name_country")) %>%
  mutate(name_region.y = name_region) %>%
  mutate(name_crop_earthstat = variable) %>%
  select("name_group_yy", "name_crop_earthstat", "name_region.y","layer_climate","score", "HDI", "landevap" )
no.dat$climate_main <- substr(no.dat$layer_climate, 1, 1)
no.dat$climate_sub <- substr(no.dat$layer_climate, 2, 2)
no.dat <- no.dat[!duplicated(no.dat)]

#Remove unusable values

exist.dat<- exist.dat[is.na(HDI) == FALSE & is.na(landevap) == FALSE]
exist.dat<- exist.dat[, -c("V1", "summass", "sumarea", "mindose", "name_chemclass", "name_country")]
exist.dat <- exist.dat[, c("name_group_yy", "name_crop_earthstat", "casrn_chemical" ,"name_region.y","wd","layer_climate","score", "HDI", "landevap" )]
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
exist.dat.cpp <- exist.dat %>%
  select(-wd)
general.scenario<- exist.dat.cpp %>%
  select(-casrn_chemical) %>%
  distinct()
#general.scenario<- general.scenario[!duplicated(general.scenario)] #Lots of replicates in scenarios, remove

#Prep dataset to fill in with ML values
#dat.out.knnlog<- matrix(data="X", nrow=length(chem.list), ncol= 5, dimnames=list(NULL, c("CASRN", "nNo", "nYes", "rsq1", "rsq2")))
dat.out.cpp <- matrix(data = "X", nrow = length(chem.list), ncol = 10, dimnames = list(NULL, c("CASRN", "nNo", "nYes","No_ActualNo", "Yes_ActualNo", "No_ActualYes",
                                                                                               "Yes_ActualYes", "No_Accuracy", "Yes_Accuracy", "CPEfficiency")))
dat.out.xgblog<- matrix(data="X", nrow=length(chem.list), ncol= 5, dimnames=list(NULL, c("CASRN", "nNo", "nYes", "rsq1", "rsq2")))
i<-3
pred_xgb_full <- data.frame()
#Building a separate model for each chem - loop through chem list and test accuracy
for(i in c(59:length(chem.list))){
  chem<- chem.list[i]
  exist.dat.chem <- exist.dat %>%
    filter(casrn_chemical == chem) #%>%
  present.dat <- exist.dat.chem
  present.dat<- present.dat %>%
    select(-casrn_chemical) %>%
    distinct()
  #present.dat<- present.dat[!duplicated(present.dat)]
  absent.dat <- general.scenario %>%
    anti_join(present.dat, by = c("score", "name_crop_earthstat", "name_region.y", 
                                  "name_group_yy", "layer_climate", "landevap", "HDI"))
  
  #Does a scenario use a chem?
  present.dat$ChemPresent<- "Yes"
  absent.dat$ChemPresent <- "No"
  absent.dat$wd<- 0
  
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
  #train.dat.cpp<- createDataPartition(machine.dat.cpp$layer_climate, p=0.7, list=FALSE)
  
  #The remaining 30% of rows are the test set - these will determine how good the model is
  #test.dat.cpp<- machine.dat.cpp[-train.dat.cpp,]
  #train.dat.cpp<- machine.dat.cpp[train.dat.cpp,]
  
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
  chemPred <- pValues.df.pred$chemPred
  #test.dat.cpp$ChemPresent <- as.character(test.dat.cpp$ChemPresent)
  #cm <- as.data.table(table(test.dat.cpp$ChemPresent, chemPred))
  
  #val.out <- c(chem,
               #nrow(absent.dat),
               #nrow(present.dat),
               #cm[V2 == "1" & chemPred == "V1"]$N,
               #cm[V2 == "1" & chemPred == "V2"]$N,
               #cm[V2 == "2" & chemPred == "V1"]$N,
               #cm[V2 == "2" & chemPred == "V2"]$N,
               #cm[V2 == "1" & chemPred == "V1"]$N/length(which(test.dat.cpp$ChemPresent == "1")),
               #cm[V2 == "2" & chemPred == "V2"]$N/length(which(test.dat.cpp$ChemPresent == "2")),
               #CPEfficiency)
  
  #dat.out.cpp[i,] <- val.out
  
  no.dat.cpp.pred <- cbind(no.dat.cpp, chemPred)
  
  #xgboost steps
  #machine.dat.xgb <- machine.dat %>%
  #filter(ChemPresent == "Yes")
  
  #only do regression model on existing chemical scenarios
  machine.dat.xgb <- present.dat
  machine.dat.xgb <- machine.dat.xgb %>%
    select(-ChemPresent) %>%
    mutate(logwd = log10(wd)) %>%
    select(-wd)
  no.dat.xgb <- no.dat.cpp.pred %>%
    filter(chemPred == "V2") %>%
    select(-ChemPresent, -chemPred)
  
  #Create a training set - this is what will be used to build the model (random 70% of machine data rows)
  #train.dat.xgb1<- createDataPartition(machine.dat.xgb$layer_climate, p=0.7, list=FALSE)
  
  #The remaining 30% of rows are the test set - these will determine how good the model is
  #test.dat.xgb<- machine.dat.xgb[-train.dat.xgb1,] 
  #train.dat.xgb<- machine.dat.xgb[train.dat.xgb1,]
  
  #xgboost to predict using caret
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
  
  newdata.xgb <- no.dat.xgb
  newdata.xgb$name_crop_earthstat[which(!(newdata.xgb$name_crop_earthstat %in% unique(machine.dat.xgb$name_crop_earthstat)))] <- NA
  newdata.xgb$name_region.y[which(!(newdata.xgb$name_region.y %in% unique(machine.dat.xgb$name_region.y)))] <- NA
  newdata.xgb$layer_climate[which(!(newdata.xgb$layer_climate %in% unique(machine.dat.xgb$layer_climate)))] <- NA
  newdata.xgb$climate_sub[which(!(newdata.xgb$climate_sub %in% unique(machine.dat.xgb$climate_sub)))] <- NA
  newdata.xgb$climate_main[which(!(newdata.xgb$climate_main %in% unique(machine.dat.xgb$climate_main)))] <- NA
  newdata.xgb$name_group_yy[which(!(newdata.xgb$name_group_yy %in% unique(machine.dat.xgb$name_group_yy)))] <- NA
  newdata <- newdata.xgb %>%
    filter(!is.na(.))
  pred.xgb <- predict(fit.model, newdata = newdata, na.action = na.pass)
  
  chem_pred_xgb_full <- cbind(newdata, pred.xgb)
  
  chem_pred_xgb_full$casrn_chemical <- chem
  
  write.csv(chem_pred_xgb_full, paste0("20220908chem_pred_xgb_full", chem, ".csv"))
  
  pred_xgb_full <- rbind(pred_xgb_full, chem_pred_xgb_full)
  
  #test.dat.xgb$name_crop_earthstat[which(!(test.dat.xgb$name_crop_earthstat %in% unique(train.dat.xgb$name_crop_earthstat)))] <- NA
  #test.dat.xgb$name_region[which(!(test.dat.xgb$name_region.y %in% unique(train.dat.xgb$name_region.y)))] <- NA
  
  # pred <- chem_pred_xgb_full$pred
  # orig <- test.dat.xgb$logwd
  # 
  # rss <- sum((pred - orig)^2)
  # tss <- sum((orig-mean(orig))^2)
  # rsq1 <- 1 - rss/tss
  # rsq2 <- cor(pred, orig)
  
  
  #knnPredict <- predict(chem_pred_knn, newdata = class_pred_test) 
  
  #confusionMatrix(knnPredict, chem_outcome_test$ChemPresent)
  
  #cm <- as.data.table(table(as.character(chem_outcome_test$ChemPresent), as.character(knnPredict)))
  
  # val.out<- c(chem,
  #             nrow(absent.dat), #Total Nos in dataset
  #             nrow(present.dat), #Total yes in dataset
  #             rsq1,
  #             rsq2)
  # 
  # dat.out.xgblog[i,]<- val.out
  
  #bootsrapping for xgboost
  #machine.dat.xgb$logwd <- as.numeric(machine.dat.xgb$logwd)
  #chem.outcome <- machine.dat.xgb %>%
    #select(logwd)
  
  # remove original variable from the data set
 # machine.dat.xgb <- machine.dat.xgb %>% 
    #select(-logwd)
  
  #str(machine.dat)
  #scale input numeric features
  #machine.dat.xgb[, c("landevap")] <- scale(machine.dat.xgb[, c("landevap")])
  #machine.dat.xgb[,c("HDI")] <- scale(machine.dat.xgb[, c("HDI")])
  #machine.dat.xgb[,c("score")] <- scale(machine.dat.xgb[, c("score")])
  #head(machine.dat)
  
  #str(machine.dat)
  #dummy code categorical features
  #name_crop_earthstat <- as.data.frame(dummy.code(machine.dat.xgb$name_crop_earthstat))
  #name_region <- as.data.frame(dummy.code(machine.dat.xgb$name_region))
  #name_cropgroup <- as.data.frame(dummy.code(machine.dat.xgb$name_cropgroup))
  #layer_climate <- as.data.frame(dummy.code(machine.dat.xgb$layer_climate))
  
  #machine.dat.xgb <- cbind(machine.dat.xgb, name_crop_earthstat,name_region, name_cropgroup,
                           #layer_climate)
  
  #machine.dat.xgb <- machine.dat.xgb %>% 
    #select(-one_of(c("name_crop_earthstat","name_region", "name_cropgroup",
                     #"layer_climate", "climate_main", "climate_sub")))
  
  #separate training & test dataset, using same separations as before
  #class_pred_train <- machine.dat.xgb[train.dat.xgb1, ]
  #label_pred_train <- as.matrix(chem.outcome[train.dat.xgb1])
  #data_pred_train <- as.matrix(cbind(class_pred_train, label_pred_train))
  
  #class_pred_test <- machine.dat.xgb[-train.dat.xgb1, ]
  #label_pred_test <- as.matrix(chem.outcome[-train.dat.xgb1])
  #data_pred_test <- as.matrix(cbind(class_pred_test, label_pred_test))
  
  #dat.xgb.rec <- 
    #recipe(logwd ~., data = data_pred_train) %>%
    #step_dummy(all_nominal())
  
  # dat.xgb.boot <- data_pred_train %>%
  #   bootstraps()
  # 
  # #dat.xgb.boot
  # 
  # dat.xgb <- 
  #   boost_tree() %>%
  #   set_mode("regression") %>%
  #   set_engine("xgboost")
  # 
  # #define a basic xgboost model
  # fit_bootstrap <- function(index){
  #   xgb_boot <- dat.xgb.boot$splits[[index]] %>%
  #     training()
  #   
  #   workflow() %>%
  #     add_recipe(dat.xgb.rec) %>%
  #     add_model(dat.xgb) %>%
  #     #boost_tree() %>%
  #     #set_mode("regression") %>%
  #     #set_engine("xgboost") %>%
  #     fit(xgb_boot) %>%
  #     write_rds(paste0("model_", index, ".rds"))
  # }
  # 
  # #fit to 25 bootstrapped datasets
  # for(j in 1:25) {
  #   fit_bootstrap(j)
  # }
  # 
  # predict_bootstrap <- function(new_data, index) {
  #   read_rds(paste0("model_", index, ".rds")) %>%
  #     predict(new_data) %>%
  #     rename(!!sym(paste0("pred_", index)) := .pred)
  # }
  # 
  # training_preds <- 
  #   seq(1,25) %>%
  #   map_dfc(~predict_bootstrap(data_pred_train, .x))
  # 
  # #training_preds
  # 
  # p <- training_preds %>%
  #   bind_cols(data_pred_train) %>%
  #   rowid_to_column() %>%
  #   pivot_longer(starts_with("pred_"),
  #                names_to = "model",
  #                values_to = ".pred") %>%
  #   group_by(rowid) %>%
  #   summarise(logwd = max(logwd),
  #             pred_mean = mean(.pred),
  #             std_dev = sd(.pred)) %>%
  #   riekelib::normal_interval(pred_mean, std_dev) %>%
  #   ggplot(aes(x = logwd,
  #              y = pred_mean)) +
  #   geom_point(alpha = 0.5) +
  #   geom_segment(aes(x = logwd,
  #                    xend = logwd,
  #                    y = ci_lower,
  #                    yend = ci_upper),
  #                alpha = 0.25)
  # 
  # #training_preds
  # p<-  seq(1,25) %>%
  #   map_dfc(~predict_bootstrap(data_pred_test, .x)) %>%
  #   bind_cols(data_pred_test) %>%
  #   rowid_to_column() %>%
  #   pivot_longer(starts_with("pred_"),
  #                names_to = "model",
  #                values_to = ".pred") %>%
  #   group_by(rowid) %>%
  #   summarise(logwd = max(logwd),
  #             pred_mean = mean(.pred),
  #             std_dev = sd(.pred)) %>%
  #   riekelib::normal_interval(pred_mean, std_dev) %>%
  #   ggplot(aes(x = logwd,
  #              y= pred_mean)) +
  #   geom_point(alpha = 0.5) +
  #   geom_segment(aes(x = logwd,
  #                    xend = logwd,
  #                    y = ci_lower,
  #                    yend = ci_upper),
  #                alpha = 0.25)
  # #ggsave(paste0("20220718pred_df_interval_testlog", chem, ".png"))
  # p <- seq(1,25) %>%
  #   map_dfc(~predict_bootstrap(data_pred_test, .x)) %>%
  #   bind_cols(data_pred_test) %>%
  #   rowid_to_column() %>%
  #   pivot_longer(starts_with("pred_"),
  #                names_to = "model",
  #                values_to = ".pred") %>%
  #   group_by(rowid) %>%
  #   summarise(logwd = max(logwd),
  #             pred_mean = mean(.pred),
  #             std_dev = sd(.pred)) %>%
  #   riekelib::normal_interval(pred_mean, std_dev)
  #write.csv(p, paste0("20220718pred_df_interval_testlog", chem, ".csv"))
  
  print(i)
}

write.csv(dat.out.cpp, "20220816dat.out.cpp_cropspecific.csv")
write.csv(dat.out.xgblog, "20220816dat.out.xgblog_cropspecific.csv")
