############################################################################

### Chem imputation interval
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
library(mosaic)
library(pracma)
#Test datasets

#exist.dat<- fread("C:/me/data gaps/20220912country_input_full_all.csv")
exist.dat<- fread("C:/me/data gaps/20221019country_input_full_all.csv")
#exist.dat.nocrop <- fread("C:/me/data gaps/20221006country_imputation_combi_existcountrynocrop.csv")
exist.dat <- exist.dat %>%
  select(-name_group_yy) %>%
  mutate(name_region.y = case_when(name_country == "CENTRAL AMERICA-CARIBBEAN" ~ "North America",
                                   TRUE ~ name_region.y)) %>%
  left_join(crop_id, by = c("name_crop_earthstat" = "CROPNAME")) %>%
  left_join(country_croparea, by = c("name_country", "name_crop_earthstat" = "variable"))
no.country.continent <- fread("C:/me/Database/harmonization/no.country.continent.csv")
#no.dat <- fread("C:/me/data gaps/20220912country_input_n_ccc.csv")
no.dat <- fread("C:/me/data gaps/20221019country_input_n_all.csv")
no.dat <- no.dat %>%
  filter(score != 0) %>%
  left_join(crop_id, by = c("variable" = "CROPNAME")) %>%
  mutate(name_group_yy = name_group_yy)%>%
  left_join(no.country.continent, by = c("name_country" = "name_country")) %>%
  left_join(country_croparea, by = c("name_country", "variable")) %>%
  mutate(name_region.y = name_region) %>%
  mutate(name_crop_earthstat = variable) %>%
  select("name_group_yy", "name_crop_earthstat", "name_region.y","layer_climate","score", "HDI", "landevap", "croparea", "sumcroparea" ) %>%
  filter(!is.na(HDI)) %>%
  filter(!is.na(landevap))
no.dat$climate_main <- substr(no.dat$layer_climate, 1, 1)
no.dat$climate_sub <- substr(no.dat$layer_climate, 2, 2)
no.dat <- no.dat[!duplicated(no.dat)]
exist.dat.nocrop <- exist.dat.nocrop %>%
  filter(score != 0) %>%
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
dat.out.cpp <- matrix(data = "X", nrow = length(chem.list), ncol = 10, dimnames = list(NULL, c("CASRN", "nNo", "nYes","No_ActualNo", "Yes_ActualNo", "No_ActualYes",
                                                                                               "Yes_ActualYes", "No_Accuracy", "Yes_Accuracy", "CPEfficiency")))
dat.out.xgblog<- matrix(data="X", nrow=length(chem.list), ncol= 5, dimnames=list(NULL, c("CASRN", "nNo", "nYes", "rsq1", "rsq2")))
i<-15
pred_xgb_full <- data.frame()
#Building a separate model for each chem - loop through chem list and test accuracy
for(i in c(1:length(chem.list))){
  chem<- chem.list[i]
  chem <- "658066-35-4"
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
  machine.dat.cpp[,1] = as.factor(machine.dat.cpp[,1])
  machine.dat.cpp[,1] = as.numeric(unlist(machine.dat.cpp[,1]))
  
  #write.csv(machine.dat.cpp, "machine.dat.cpp_658066-35-4.csv")
  
  no.dat$ChemPresent = "null"
  no.dat.cpp <- no.dat %>%
    select(ChemPresent, name_crop_earthstat, name_region.y, name_group_yy, HDI, landevap, layer_climate, climate_main, climate_sub, score, croparea, sumcroparea)
  #write.csv(no.dat.cpp, "no.dat.cpp_658066-35-4.csv")
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
  
  #cpp bootstrapping
  #n <- nrow(machine.dat.cpp)
  #nbootstraps <- as.integer(sqrt(n))
  #bootstrap_preds_no <- vector(length = nbootstraps)
  #bootstrap_preds_yes <- vector(length = nbootstraps)
  #bootstrap_preds_no <- vector(length = 10)#remember to change the number to nbootstraps
  #bootstrap_preds_yes <- vector(length = 10)#remember to change the number to nbootstraps
  bootstrap_preds_yes <- matrix(data = "X", nrow = nrow(no.dat.cpp), ncol = 10)
  bootstrap_preds_no <- matrix(data = "X", nrow = nrow(no.dat.cpp), ncol = 10)
  #bootstrapping
  #cpp_boot <- matrix(data="X", nrow=nrow(machine.dat.cpp), ncol= 4,  dimnames=list(NULL, c("yesbootlower","yesbootupper", "nobootlower","nobootupper")))
  #for (i in c(1:10)){
  for(b in c(1:10)){
    #train_idxs = sample(seq(n),size = n, replace = TRUE)
    originalData_r <- machine.dat.cpp
    #originalData_r <- originalData_r %>%
      #select(-orig.id)
    pValues = ICPClassification(originalData_r, no.dat.cpp, ratioTrain = 0.7, method = "rf", nrTrees = 100)
    for(i in c(1:nrow(no.dat.cpp))){
      bootstrap_preds_no[i,b] = pValues[i,1]
      bootstrap_preds_yes[i,b] = pValues[i,2]
      #print(i)
    }
    print(b)
  }
  
  # bootstrap_preds_no.df <- as.data.frame(bootstrap_preds_no)
  # bootstrap_preds_no.df10 <- bootstrap_preds_no.df[1:20,]
  # bootstrap_preds_no.df <- bootstrap_preds_no.df %>%
  #   mutate_if(is.character, as.numeric) #%>%
  #   rowwise() %>%
  #   mutate(median = median(c_across(where(is.numeric)), na.rm = TRUE))
  # percentile_boot_no <- bootstrap_preds_no.df %>%
  #   bind_cols(as.data.frame(t(apply(., 1, quantile, c(0.05, 0.95)))))
  # percentile_boot_no <- percentile_boot_no %>%
  #   mutate(no0.05 = `5%`, no0.95 = `95%`) %>%
  #   dplyr::select(no0.05, no0.95)
  # bootstrap_preds_no_median <- bootstrap_preds_no.df %>%
  #   rowwise() %>%
  #   mutate(median = median(c_across(where(is.numeric)), na.rm = TRUE))
  # percentile_median_no <- cbind(percentile_boot_no, bootstrap_preds_no_median$median)
  
  
  
  bootstrap_preds_yes.df <- as.data.frame(bootstrap_preds_yes)
  bootstrap_preds_yes.df <- bootstrap_preds_yes.df %>%
    mutate_if(is.character, as.numeric)
  percentile_boot_yes10 <- bootstrap_preds_yes.df %>%
    bind_cols(as.data.frame(t(apply(., 1, quantile, c(0.025, 0.975)))))
  percentile_boot_yes <- percentile_boot_yes10 %>%
    mutate(yes0.025 = `2.5%`, yes0.975 = `97.5%`) %>%
    select(yes0.025, yes0.975)
  bootstrap_preds_yes_median <- bootstrap_preds_yes.df %>%
    rowwise() %>%
    mutate(median = median(c_across(where(is.numeric)), na.rm = TRUE))
  percentile_median_yes <- cbind(no.dat.cpp, percentile_boot_yes, bootstrap_preds_yes_median) %>%
    filter(median > 0.05)
  write.csv(percentile_median_yes, "20221020percentile_median_yes.csv")
  
  
  percentile_median_yes %>%
    arrange(`bootstrap_preds_yes_median$median`) %>%
    mutate(id = rownames(.)) %>%
    ggplot(aes(x = reorder(id, `bootstrap_preds_yes_median$median`))) +
    geom_point(aes(y = `bootstrap_preds_yes_median$median`)) +
    geom_errorbar(aes(ymin=yes0.025,
                      ymax=yes0.975),
                  width=0.01, color = "grey20", alpha = 0.5) +
    #geom_errorbar(aes(ymin=CIlow,
    #ymax=CIhigh),
    #width=0.01, color = "red", alpha = 0.2) +
    #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
    #labels = trans_format("log10", math_format(10^.x)))+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    geom_hline(yintercept=0.05, linetype="dashed", color = "red")
  ggsave(paste0("20221019cpp.yes.interval",chem,".jpeg"), width = 20, height = 20, units = "cm")
  #percentile_yesno <- cbind(percentile_boot_no, percentile_boot_yes)
  
  #no.dat.cpp.pred <- cbind(no.dat.cpp, chemPred)
  
  #no.dat.cpp.pred.interval <- cbind(no.dat.cpp.pred, pValues.df,percentile_yesno)
  #no.dat.cpp.pred.interval <- no.dat.cpp.pred.interval %>%
    #filter(chemPred == "V2") %>%
    #select(-ChemPresent, -chemPred)
  
  #xgboost steps
  #machine.dat.xgb <- machine.dat %>%
  #filter(ChemPresent == "Yes")
  
  #only do regression model on existing chemical scenarios
  # machine.dat.xgb <- present.dat
  # machine.dat.xgb <- machine.dat.xgb %>%
  #   select(-ChemPresent) %>%
  #   mutate(logwd = log10(wd)) %>%
  #   select(-wd)
  # no.dat.xgb <- no.dat.cpp.pred %>%
  #   filter(chemPred == "V2") %>%
  #   select(-ChemPresent, -chemPred)
  
  #Create a training set - this is what will be used to build the model (random 70% of machine data rows)
  #train.dat.xgb1<- createDataPartition(machine.dat.xgb$layer_climate, p=0.7, list=FALSE)
  
  #The remaining 30% of rows are the test set - these will determine how good the model is
  #test.dat.xgb<- machine.dat.xgb[-train.dat.xgb1,] 
  #train.dat.xgb<- machine.dat.xgb[train.dat.xgb1,]
  
  # #xgboost to predict using caret
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
  # newdata.xgb <- no.dat.xgb
  # newdata.xgb$name_crop_earthstat[which(!(newdata.xgb$name_crop_earthstat %in% unique(machine.dat.xgb$name_crop_earthstat)))] <- NA
  # newdata.xgb$name_region.y[which(!(newdata.xgb$name_region.y %in% unique(machine.dat.xgb$name_region.y)))] <- NA
  # newdata.xgb$layer_climate[which(!(newdata.xgb$layer_climate %in% unique(machine.dat.xgb$layer_climate)))] <- NA
  # newdata.xgb$climate_sub[which(!(newdata.xgb$climate_sub %in% unique(machine.dat.xgb$climate_sub)))] <- NA
  # newdata.xgb$climate_main[which(!(newdata.xgb$climate_main %in% unique(machine.dat.xgb$climate_main)))] <- NA
  # newdata.xgb$name_group_yy[which(!(newdata.xgb$name_group_yy %in% unique(machine.dat.xgb$name_group_yy)))] <- NA
  # newdata <- newdata.xgb %>%
  #   filter(!is.na(.))
  # newdata <- newdata %>%
  #   filter(!is.na(layer_climate)) %>%
  #   filter(!is.na(name_crop_earthstat))
  # pred.xgb <- predict(fit.model, newdata = newdata, na.action = na.pass)
  # 
  # chem_pred_xgb_full <- cbind(newdata, pred.xgb)
  # 
  # chem_pred_xgb_full$casrn_chemical <- chem
  # 
  # write.csv(chem_pred_xgb_full, paste0("20220906chem_pred_xgb_full", chem, ".csv"))
  # 
  # pred_xgb_full <- rbind(pred_xgb_full, chem_pred_xgb_full)
  # 
  # 
  # #bootsrapping for xgboost
  # n <- nrow(machine.dat.xgb)
  # nbootstraps <- as.integer(sqrt(n))
  # bootstrap_preds <- vector(length = nbootstraps)
  # bootstrap_preds <- vector(length = 5) #change number
  # val_residuals <- vector()
  # preds.xgb.boot.full <- no.dat.xgb
  # preds.xgb.boot <- matrix(data = "X", nrow = nrow(no.dat.xgb), ncol = 5)
  # xgb.quantile.boot <- matrix(data = "X", nrow = nrow(no.dat.xgb), ncol = 4)
  # #xgb_boot <- matrix(data="X", nrow=nrow(bbch.n), ncol= 12,  dimnames=list(NULL, c("predyes","yesbootlower","yesbootupper", "predno","nobootlower","nobootupper","generalisationyes", "yeslower", "yesupper",  "generalisationno","nolower", "noupper")))
  # #for (i in c(31:100)){
  # write.csv(no.dat.xgb, "no.dat.xgb.csv")
  # no.dat.xgb <- read.csv("no.dat.xgb.csv")
  # no.dat.xgb <- no.dat.xgb %>%
  #   select(-X)
  # preds.xgb.boot.full <- no.dat.xgb
  # for(b in c(1:20)){
  #     train_idxs <- sample(seq(n),size = n, replace = TRUE)
  #     fit.model <- train(logwd ~., data = machine.dat.xgb[train_idxs,], method = "xgbTree", metric = "RMSE", trControl = control)
  #     newdata.xgb <- no.dat.xgb
  #     newdata.xgb$name_crop_earthstat[which(!(newdata.xgb$name_crop_earthstat %in% unique(machine.dat.xgb[train_idxs,]$name_crop_earthstat)))] <- NA
  #     newdata.xgb$name_region.y[which(!(newdata.xgb$name_region.y %in% unique(machine.dat.xgb[train_idxs,]$name_region.y)))] <- NA
  #     newdata.xgb$layer_climate[which(!(newdata.xgb$layer_climate %in% unique(machine.dat.xgb[train_idxs,]$layer_climate)))] <- NA
  #     newdata.xgb$climate_sub[which(!(newdata.xgb$climate_sub %in% unique(machine.dat.xgb[train_idxs,]$climate_sub)))] <- NA
  #     newdata.xgb$climate_main[which(!(newdata.xgb$climate_main %in% unique(machine.dat.xgb[train_idxs,]$climate_main)))] <- NA
  #     newdata.xgb$name_group_yy[which(!(newdata.xgb$name_group_yy %in% unique(machine.dat.xgb[train_idxs,]$name_group_yy)))] <- NA
  #     newdata <- newdata.xgb %>%
  #       filter(!is.na(.))
  #     newdata.val <- machine.dat.xgb[-train_idxs,] 
  #     newdata.val$name_crop_earthstat[which(!(newdata.val$name_crop_earthstat %in% unique(machine.dat.xgb[train_idxs,]$name_crop_earthstat)))] <- NA
  #     newdata.val$name_region.y[which(!(newdata.val$name_region.y %in% unique(machine.dat.xgb[train_idxs,]$name_region.y)))] <- NA
  #     newdata.val$layer_climate[which(!(newdata.val$layer_climate %in% unique(machine.dat.xgb[train_idxs,]$layer_climate)))] <- NA
  #     newdata.val$climate_sub[which(!(newdata.val$climate_sub %in% unique(machine.dat.xgb[train_idxs,]$climate_sub)))] <- NA
  #     newdata.val$climate_main[which(!(newdata.val$climate_main %in% unique(machine.dat.xgb[train_idxs,]$climate_main)))] <- NA
  #     newdata.val$name_group_yy[which(!(newdata.val$name_group_yy %in% unique(machine.dat.xgb[train_idxs,]$name_group_yy)))] <- NA
  #     newdata.val <- newdata.val %>%
  #       filter(!is.na(.))
  #     preds <- predict(fit.model, newdata = newdata.val, na.action = na.pass)
  #     val_residuals <- append(val_residuals,(preds - newdata.val$logwd))
  #     preds.xgb <- predict(fit.model, newdata = newdata, na.action = na.pass)
  #     preds.xgb.boot <- cbind(newdata, preds.xgb)
  #     varname  <- paste0("preds.xgb", b)
  #     preds.xgb.boot[[varname]] <- preds.xgb.boot$preds.xgb
  #     write.csv(preds.xgb.boot, "preds.xgb.boot.csv")
  #     preds.xgb.boot <- read.csv("preds.xgb.boot.csv")
  #     preds.xgb.boot <- preds.xgb.boot %>%
  #       select(-X) %>%
  #       select(-preds.xgb)
  #     preds.xgb.boot.full <- preds.xgb.boot.full %>%
  #       left_join(preds.xgb.boot, by = c("name_crop_earthstat", "name_region.y", "name_group_yy", "HDI", "landevap", "layer_climate", "climate_main", "climate_sub", "score"))
  #     print(b)
  #  }
  # preds.xgb.boot.full.rowmean <- preds.xgb.boot.full %>%
  #   mutate(mean = rowMeans(.[10:29], na.rm = TRUE))%>%
  #   mutate(across(preds.xgb1:preds.xgb5, ~ . -mean))
  # #preds.xgb.mean = preds.xgb - mean(preds.xgb)
  # val_residuals <- cbind(val_residuals)
  # 
  # fit.model <- train(logwd ~., data = machine.dat.xgb, method = "xgbTree", metric = "RMSE", trControl = control)
  # preds <- predict(fit.model, newdata = machine.dat.xgb, na.action = na.pass)
  # train_residuals <- machine.dat.xgb$logwd - preds
  # 
  # val_residuals <- quantile(val_residuals, seq(0,1,by = 0.01))
  # train_residuals <- quantile(train_residuals, seq(0,1,by = 0.01))
  # 
  # alpha <- 0.05
  # qs = c(alpha/2, 1-alpha/2)  
  # 
  # no_information_error <- mean(abs(randperm(machine.dat.xgb$logwd) - randperm(preds)))
  # generalisation <- abs(mean(val_residuals) - mean(train_residuals))
  #   
  # no_information_val <- abs(no_information_error - train_residuals)
  # relative_overfitting_rate <- mean(generalisation/no_information_val)
  #   
  # weight <- 0.632 / (1 - 0.368 * relative_overfitting_rate)
  # residuals <- (1-weight) * train_residuals + weight * val_residuals
  # 
  # for(a in c(1:nrow(preds.xgb.boot.full))){
  #   C = expand.grid(a= preds.xgb.boot.full[a, 10:14], b = residuals)$a + expand.grid(a= preds.xgb.boot.full[a, 10:14], b = residuals)$b
  #   percentile = quantile(C, qs, na.rm = TRUE)
  #   xgb.quantile.boot[a, 1] <- percentile[1] 
  #   xgb.quantile.boot[a,2] <- percentile[2]
  #   xgb.quantile.boot[a,3] <- quantile(preds.xgb.boot.full[a, 10:14], 0.025, na.rm = TRUE)
  #   xgb.quantile.boot[a,4] <- quantile(preds.xgb.boot.full[a, 10:14], 0.975, na.rm = TRUE)
  #   print(a)
  # }
  

  
  print(i)
}


fit.model <- train(logwd ~., data = machine.dat.xgb, method = "xgbTree", metric = "RMSE", trControl = control)
newdata.xgb <- no.dat.xgb
newdata.xgb$name_crop_earthstat[which(!(newdata.xgb$name_crop_earthstat %in% unique(machine.dat.xgb$name_crop_earthstat)))] <- NA
newdata.xgb$name_region.y[which(!(newdata.xgb$name_region.y %in% unique(machine.dat.xgb$name_region.y)))] <- NA
newdata.xgb$layer_climate[which(!(newdata.xgb$layer_climate %in% unique(machine.dat.xgb$layer_climate)))] <- NA
newdata.xgb$climate_sub[which(!(newdata.xgb$climate_sub %in% unique(machine.dat.xgb$climate_sub)))] <- NA
newdata.xgb$climate_main[which(!(newdata.xgb$climate_main %in% unique(machine.dat.xgb$climate_main)))] <- NA
newdata.xgb$name_group_yy[which(!(newdata.xgb$name_group_yy %in% unique(machine.dat.xgb$name_group_yy)))] <- NA
newdata <- newdata.xgb %>%
  filter(!is.na(.))
preds <- predict(fit.model, newdata = newdata, na.action = na.pass)
newdata.preds <- cbind(newdata, preds)
preds.xgb.full <- no.dat.xgb %>%
  left_join(newdata.preds, by = c("name_crop_earthstat", "name_region.y", "name_group_yy", "HDI", "landevap", "layer_climate", "climate_main", "climate_sub", "score"))

preds.xgb.full.interval %>%
  rowwise() %>%
  filter(`2` < `4`)

preds.xgb.full.interval <- cbind(preds.xgb.full, xgb.quantile.boot)
write.csv(preds.xgb.full.interval, "20220906preds.xgb.full.interval.csv")

preds.cppxgb.full.interval <- no.dat.cpp.pred.interval %>%
  full_join(preds.xgb.full.interval, by = c("name_crop_earthstat", "name_region.y", "name_group_yy", "HDI", "landevap", "layer_climate", "climate_main", "climate_sub", "score"))

write.csv(preds.cppxgb.full.interval, "20220907preds.cppxgb.full.interval.csv")


write.csv(dat.out.cpp, "20220816dat.out.cpp_cropspecific.csv")
write.csv(dat.out.xgblog, "20220816dat.out.xgblog_cropspecific.csv")


preds.xgb.full.plot <- preds.xgb.boot.full[9:14, 10:29]
preds.xgb.full.plot <- data.frame(t(preds.xgb.full.plot))
preds.xgb.full.plot %>%
  mutate(dose = 10^(X13)) %>%
  ggplot(aes(x = dose)) +
  geom_histogram() +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))

preds.xgb.full.interval %>%
  mutate(preds = 10^as.numeric(preds)) %>%
  mutate(CIlow_adj = 10^as.numeric(`1`)) %>%
  mutate(CIhigh_adj = 10^as.numeric(`2`)) %>%
  mutate(CIlow = 10^as.numeric(`3`)) %>%
  mutate(CIhigh = 10^as.numeric(`4`)) %>%
  filter(!is.na(preds)) %>%
  arrange(preds) %>%
  mutate(id = rownames(.)) %>%
  ggplot(aes(x = reorder(id, preds))) +
  geom_point(aes(y = preds)) +
  geom_errorbar(aes(ymin=CIlow_adj,
                    ymax=CIhigh_adj),
                width=0.01, color = "grey20", alpha = 0.5) +
  geom_errorbar(aes(ymin=CIlow,
                    ymax=CIhigh),
                width=0.01, color = "red", alpha = 0.2) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("20220907preds.xgb.full.interval2.jpeg")

percentile_median_yes %>%
  arrange(`bootstrap_preds_yes_median$median`) %>%
  mutate(id = rownames(.)) %>%
  ggplot(aes(x = reorder(id, `bootstrap_preds_yes_median$median`))) +
  geom_point(aes(y = `bootstrap_preds_yes_median$median`)) +
  geom_errorbar(aes(ymin=yes0.05,
                    ymax=yes0.95),
                width=0.01, color = "grey20", alpha = 0.5) +
  #geom_errorbar(aes(ymin=CIlow,
                    #ymax=CIhigh),
                #width=0.01, color = "red", alpha = 0.2) +
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                #labels = trans_format("log10", math_format(10^.x)))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_hline(yintercept=0.05, linetype="dashed", color = "red")
ggsave("20221012cpp.yes.interval.jpeg", width = 20, height = 20, units = "cm")
