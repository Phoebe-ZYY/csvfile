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

setwd("C:/Users/markos/Documents/DataDownloads/BBCH_Sep_Files")
#Existed: complete info including chem use for crop-country combos in EU 
bbch.dat<- fread("bbch_input.csv")


bbch.dat <- fread("C:/me/data gaps/bbch_input_full.csv")
bbch.impute <- fread("C:/me/data gaps/bbch_impute_specific.csv")

head(bbch.dat)
summary(bbch.dat)
table(bbch.dat$name_cropgroup.y)
table(bbch.dat$name_crop_earthstat)
table(bbch.dat$name_ai.x)
table(bbch.dat$name_country)
table(bbch.dat$class1)
table(bbch.dat$layer_climate)
table(bbch.dat$name_indication)
quantile(bbch.dat$gpw_2015) #continuous


test<- bbch.dat[is.na(HDI_2015) == FALSE & is.na(GDP_2015) == FALSE & is.na(landevap_2017) == FALSE & 
                   is.na(tmin_all_mean) == FALSE & is.na(tmax_all_mean) == FALSE & is.na(precip_all_mean) == FALSE & 
                   is.na(ocs) == FALSE & is.na(soc_15_proj) == FALSE & is.na(soc_5_proj) == FALSE]
test<- bbch.dat[is.na(HDI) == FALSE & is.na(landevap) == FALSE]
impute<- bbch.impute[is.na(HDI) == FALSE & is.na(landevap) == FALSE]

nrow(test)/nrow(bbch.dat) #work with this 84% for now
bbch.dat<- test

nrow(impute)/nrow(bbch.impute) #work with this 84% for now
bbch.impute<- impute


quantile(bbch.dat$tmin_all_mean)
quantile(bbch.dat$precip_all_mean)
quantile(bbch.dat$HDI_2015)
quantile(bbch.dat$GDP_2015)
quantile(bbch.dat$landevap_2017)
quantile(bbch.dat$ocs) #continuous
quantile(bbch.dat$soc_15_proj) #continuous
quantile(bbch.dat$soc_5_proj)  #continuous

#Try developing separate models for each bbch

#Check data types again
summary(bbch.dat)

#Remove unusable values
length(unique(bbch.dat$name_ai)) == length(unique(bbch.dat$casrn_chemical)) #Only use 1 of these columns

bbch.dat<- bbch.dat[, -c("V1", "name_ai", "dose", "casrn_chemical")]
bbch.impute<- bbch.impute[, -c("V1", "name_ai", "dose")]
bbch.impute.f <- bbch.impute[,-c("casrn_chemical")]
bbch.impute.f <- bbch.impute.f[!duplicated(bbch.impute.f)]

#Check on different input data present
table(bbch.dat$bbch_specific)

#Make sure no duplicates
bbch.dat<- bbch.dat[!duplicated(bbch.dat)]

#Now build ML model - still many options for each, so try without large classes for now
bbch.res<- bbch.dat
bbch.dat<- bbch.dat[, list(name_cropgroup.y, name_indication, bbch_specific, layer_climate,
                        HDI_2015, GDP_2015, landevap_2017, gpw_2015, tmin_all_mean, tmax_all_mean, precip_all_mean,
                        ocs, soc_15_proj, soc_5_proj)]

bbch.list<- unique(bbch.dat$bbch_specific)
general.scenario<- bbch.dat[, -c("bbch_specific")]

general.scenario<- general.scenario[!duplicated(general.scenario)] #Lots of replicates in scenarios

dat.out<- matrix(data="X", nrow=length(bbch.list), ncol= 9, dimnames=list(NULL, c("bbch_specific", "nNo", "nYes", "No_ActualNo", "Yes_ActualNo",  
                                                                                  "No_ActualYes", "Yes_ActualYes", "No_Accuracy", "Yes_Accuracy")))

for(i in c(1:length(bbch.list))){
  bbch<- bbch.list[i]
  present.dat<- bbch.dat[bbch_specific %in% bbch] #ID scenarios that use that chem
  present.dat<- present.dat[, -c("bbch_specific")]
  present.dat<- present.dat[!duplicated(present.dat)]
  x<- row.match(present.dat, general.scenario) #ID scenarios that do not use that chem
  absent.dat<- general.scenario[!na.omit(x)]
  
  test<- row.match(present.dat, absent.dat)
  if(length(which(is.na(test) == FALSE))) stop ('Multi-value match')
  
  #Does a scenario use a chem?
  present.dat$ChemPresent<- "Yes"
  absent.dat$ChemPresent<- "No"
  
  machine.dat<- as.data.table(rbind(present.dat, absent.dat))
  
  machine.dat<- transform(
    machine.dat,
    #name_country=as.factor(as.character(name_country)),
    #name_crop_earthstat=as.factor(as.character(name_crop_earthstat)),
    name_cropgroup.y=as.factor(as.character(name_cropgroup.y)),
    #casrn_chemical=as.factor(as.character(casrn_chemical)),
    #class1=as.factor(as.character(class1)),
    name_indication=as.factor(as.character(name_indication)),
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
    soc_5_proj=as.numeric(as.character(soc_5_proj)),
    ChemPresent=as.factor(ChemPresent)
  )
  
  train.dat<- createDataPartition(machine.dat$ChemPresent, p=0.7, list=FALSE)
  test.dat<- machine.dat[-train.dat,] 
  train.dat<- machine.dat[train.dat,] 
  
  rf <- randomForest(
    ChemPresent ~ .,
    data=train.dat
  )
  
  pred <- predict(rf, newdata=test.dat[,-"ChemPresent"])
  cm <- as.data.table(table(test.dat$ChemPresent, pred))
  
  val.out<- c(bbch,
              nrow(absent.dat), #Total Nos in dataset
              nrow(present.dat), #Total yes in dataset
              cm[V2=="No" & pred == "No"]$N, #No actual no
              cm[V2=="No" & pred == "Yes"]$N, #Yes actual no
              cm[V2=="Yes" & pred == "No"]$N, #No actual yes
              cm[V2=="Yes" & pred == "Yes"]$N, #Yes actual yes
              cm[V2=="No" & pred == "No"]$N/length(which(test.dat$ChemPresent == "No")), #No accuracy
              cm[V2=="Yes" & pred == "Yes"]$N/length(which(test.dat$ChemPresent == "Yes")) #Yes accuracy
              ) 
  
  dat.out[i,]<- val.out
  print(i)
}
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

setwd("C:/Users/markos/Documents/DataDownloads/BBCH_Sep_Files")
#save.image(file="210915-Preliminary_ML_Analysis.RData")
save.image(file="220112-Sep_BBCH_Analysis.RData")
#load("220112-Sep_BBCH_Analysis.RData")

bbch<- "1092"
present.dat<- bbch.dat[bbch_specific %in% bbch] #ID scenarios that use that chem
present.dat<- present.dat[, -c("bbch_specific")]
present.dat<- present.dat[!duplicated(present.dat)]
x<- row.match(present.dat, general.scenario) #ID scenarios that do not use that chem
absent.dat<- general.scenario[!na.omit(x)]

test<- row.match(present.dat, absent.dat)
if(length(which(is.na(test) == FALSE))) stop ('Multi-value match')

#Does a scenario use a chem?
present.dat$ChemPresent<- "Yes"
absent.dat$ChemPresent<- "No"

machine.dat<- as.data.table(rbind(present.dat, absent.dat))

machine.dat<- transform(
  machine.dat,
  name_country=as.factor(as.character(name_country)),
  name_crop_earthstat=as.factor(as.character(name_crop_earthstat)),
  name_cropgroup=as.factor(as.character(name_cropgroup)),
  #casrn_chemical=as.factor(as.character(casrn_chemical)),
  #class1=as.factor(as.character(class1)),
  name_indication=as.factor(as.character(name_indication)),
  layer_climate=as.factor(as.character(layer_climate)),
  HDI=as.numeric(as.character(HDI)),
  #GDP_2015=as.numeric(as.character(GDP_2015)),
  landevap=as.numeric(as.character(landevap)),
  name_chemclass1 = as.character(name_chemclass1),
  name_moagroup = as.character(name_moagroup),
  #gpw_2015=as.numeric(as.character(gpw_2015)),
  #tmin_all_mean=as.numeric(as.character(tmin_all_mean)),
  #tmax_all_mean=as.numeric(as.character(tmax_all_mean)),
  #precip_all_mean=as.numeric(as.character(precip_all_mean)),
  #ocs=as.numeric(as.character(ocs)),
  #soc_15_proj=as.numeric(as.character(soc_15_proj)),
  #soc_5_proj=as.numeric(as.character(soc_5_proj)),
  ChemPresent=as.factor(ChemPresent)
)
write.csv(machine.dat, "C:/me/data gaps/bbch_1092.csv")
write.csv(bbch.impute.f, "C:/me/data gaps/bbch_impute_f.csv")
