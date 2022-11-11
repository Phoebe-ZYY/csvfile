############################################################################

### Chemical Separate Models knn
###

############################################################################

rm(list = ls())

library(data.table)
library(prodlim)
library(randomForest)
library(caret)
library(class)
library(dplyr)
library(e1071)
library(FNN) 
library(gmodels) 
library(psych)
#Sys.setenv(LANG = "en")

#Test datasets

exist.dat <- fread("C:/me/data gaps/20220609country_input_full.csv")
exist.dat$climate_main <- substr(exist.dat$layer_climate, 1, 1)
exist.dat$climate_sub <- substr(exist.dat$layer_climate, 2, 2)

exist.dat<- exist.dat[, -c("V1", "summass", "sumarea", "mindose", "name_chemclass", "name_country")]

exist.dat <- exist.dat[is.na(HDI) == FALSE & is.na(landevap) == FALSE]

chem.list<- unique(exist.dat$casrn_chemical)

#nrow(test)/nrow(bbch.dat) #work with this 84% for now
#bbch.dat<- test




#Remove unusable values
#length(unique(bbch.dat$name_ai.x)) == length(unique(bbch.dat$casrn_chemical)) #Only use 1 of these columns

#bbch.dat<- bbch.dat[, -c("V1", "name_ai.x")]

#Check on different input data present
#table(bbch.dat$bbch_specific)

#Prep set of general scenarios that match chemical-specific scenarios
general.scenario<- exist.dat[, -c("casrn_chemical")]
general.scenario<- general.scenario[!duplicated(general.scenario)] #Lots of replicates in scenarios, remove

dat.out<- matrix(data="X", nrow=length(chem.list), ncol= 2, dimnames=list(NULL, c("CASRN", "rsq")))
i <- 1
#for(i in c(1:length(chem.list))){
for(i in c(11:50)){  
  chem<- chem.list[i]
  present.dat<- exist.dat[casrn_chemical %in% chem] #ID scenarios that use that chem
  present.dat<- present.dat[, -c("casrn_chemical")]
  present.dat<- present.dat[!duplicated(present.dat)]
  #x<- row.match(present.dat, general.scenario) #ID scenarios that do not use that chem
  #absent.dat<- general.scenario[!na.omit(x)]
  absent.dat <- general.scenario %>%
    anti_join(present.dat, by = c("score", "name_crop_earthstat", "name_region",
                                  "name_cropgroup", "wd", "layer_climate", 
                                  "landevap", "HDI"))
  
  #Do any of the present rows match absent rows? If yes, there's an error
  #test<- row.match(present.dat, absent.dat)
  #if(length(which(is.na(test) == FALSE)) != 0) stop ('Multi-value match')
  
  #Does a scenario have that chemical?
  present.dat$ChemPresent<- "Yes"
  absent.dat$ChemPresent<- "No"
  absent.dat$wd <- 0
  
  machine.dat<- as.data.table(rbind(present.dat, absent.dat))
  
  # put outcome in its own object
  #machine.dat$ChemPresent <- as.factor(machine.dat$ChemPresent)
  #chem.outcome <- machine.dat %>% 
    #select(ChemPresent)
  
  machine.dat$wd <- as.numeric(machine.dat$wd)
  chem.outcome <- machine.dat %>%
    select(wd)
  
  # remove original variable from the data set
  machine.dat <- machine.dat %>% 
    select(-wd)
  
  #str(machine.dat)
  machine.dat[, c("landevap")] <- scale(machine.dat[, c("landevap")])
  machine.dat[,c("HDI")] <- scale(machine.dat[, c("HDI")])
  machine.dat[,c("score")] <- scale(machine.dat[, c("score")])
  #head(machine.dat)
  
  #str(machine.dat)
  
  #name_crop_earthstat <- as.data.frame(dummy.code(machine.dat$name_crop_earthstat))
  name_region <- as.data.frame(dummy.code(machine.dat$name_region))
  name_cropgroup <- as.data.frame(dummy.code(machine.dat$name_cropgroup))
  layer_climate <- as.data.frame(dummy.code(machine.dat$layer_climate))
  ChemPresent <- as.data.frame(dummy.code(machine.dat$ChemPresent))
  #climate_main <- as.data.frame(dummy.code(machine.dat$climate_main))
  #climate_sub <- as.data.frame(dummy.code(machine.dat$name_climate_sub))
  
  #without crop_earthstat
  #machine.dat <- cbind(machine.dat, name_crop_earthstat, name_region, name_cropgroup,
                       #layer_climate, climate_main)
  machine.dat <- cbind(machine.dat, name_region, name_cropgroup,
                       layer_climate, ChemPresent)
  
  machine.dat <- machine.dat %>% 
    select(-one_of(c("name_crop_earthstat", "name_region", "name_cropgroup",
                     "layer_climate", "climate_main", "climate_sub", "ChemPresent")))
  #head(machine.dat)
  
  set.seed(1234) # set the seed to make the partition reproducible
  
  # 0% of the sample size
  smp_size <- floor(0.7 * nrow(machine.dat))
  
  train_ind <- sample(seq_len(nrow(machine.dat)), size = smp_size)
  
  # creating test and training sets that contain all of the predictors
  class_pred_train <- machine.dat[train_ind, ]
  class_pred_test <- machine.dat[-train_ind, ]
  
  chem_outcome_train <- chem.outcome[train_ind, ]
  chem_outcome_test <- chem.outcome[-train_ind, ]
  
  #knn
  chem_pred_knn <- FNN::knn.reg(train = class_pred_train, 
                       test = class_pred_test, y = chem_outcome_train$wd, k=9)
  
  # put "mjob_outcome_test" in a data frame
  #chem_outcome_test <- data.frame(chem_outcome_test)
  
  # merge "mjob_pred_knn" and "mjob_outcome_test" 
  #class_comparison <- data.frame(chem_pred_knn, chem_outcome_test)
  
  # specify column names for "class_comparison"
  #names(class_comparison) <- c("Predictedchem", "Observedchem")
  
  # inspect "class_comparison" 
  #head(class_comparison)
  
  # create table examining model accuracy
  #CrossTable(x = class_comparison$Observedchem, y = class_comparison$Predictedchem, 
             #prop.chisq=FALSE, prop.c = FALSE, prop.r = FALSE, prop.t = FALSE)
  
  #caret
  #chem_pred_caret <- train(class_pred_train, chem_outcome_train$wd, 
                           #method = "knn", preProcess = c("center","scale"))
  #chem_pred_caret
  #plot(chem_pred_caret)
  
  Predict_knn <- as.data.frame(chem_pred_knn$pred)
  chem_pred_knn_full <- cbind(class_pred_test, chem_outcome_test,Predict_knn)
  
  chem_pred_knn_full %>%
    ggplot(aes(x = wd, y = `chem_pred_knn$pred`)) +
    geom_point() +
    geom_abline()
  ggsave(paste0("20220622pred_df_knn_chempresent", chem, ".png"))
  
  write.csv(chem_pred_knn_full, paste0("20220622pred_df_knn_chempresent", chem, ".csv"))
  
  pred <- chem_pred_knn_full$`chem_pred_knn$pred`
  orig <- chem_pred_knn_full$wd
  
  rsq <- cor(orig, pred) ^ 2
  
  #knnPredict <- predict(chem_pred_knn, newdata = class_pred_test) 
  
  #confusionMatrix(knnPredict, chem_outcome_test$ChemPresent)
  
  #cm <- as.data.table(table(as.character(chem_outcome_test$ChemPresent), as.character(knnPredict)))
  
  val.out<- c(chem,
              rsq
  ) 
  
  dat.out[i,]<- val.out
  print(i)
}

dat.output <- as.data.table(dat.out) 
fwrite(dat.output, file = "20220622-chemical_seperate_knn50.csv")



fwrite(var.dat.out, file="220416-Full_Variables_Sep_BBCH_RFE300.csv")
fwrite(cm.dat.out, file="220416-Full_CM_Sep_BBCH_RFE300.csv")
fwrite(imp.dat.out, file="220416-Full_Important_Variables_Sep_BBCH_RFE300.csv")
fwrite(res.dat.out, file="220416-Full_Results_Sep_BBCH_RFE300.csv")
