############################################################################

### BBCH Separate Models RFE
###
### Author: Marissa Kosnik

############################################################################

rm(list = ls())

library(data.table)
library(prodlim)
library(caret)
library(randomForest)

Sys.setenv(LANG = "en")

#Test datasets

setwd("C:/Users/markos/Documents/DataDownloads/BBCH_Sep_Files")
#Existed: complete info including chem use for crop-country combos in EU 
bbch.dat<- fread("bbchspecific_input.csv")

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

nrow(test)/nrow(bbch.dat) #work with this 84% for now
bbch.dat<- test

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
length(unique(bbch.dat$name_ai.x)) == length(unique(bbch.dat$casrn_chemical)) #Only use 1 of these columns

bbch.dat<- bbch.dat[, -c("V1", "name_ai.x")]

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

#Now do rfe for each split and see what happens

for(i in c(1:length(bbch.list))){
  bbch<- bbch.list[i]
  present.dat<- bbch.dat[bbch_specific %in% bbch] #ID scenarios that use that chem
  present.dat<- present.dat[, -c("bbch_specific")]
  present.dat<- present.dat[!duplicated(present.dat)]
  x<- row.match(present.dat, general.scenario) #ID scenarios that do not use that chem
  absent.dat<- general.scenario[!na.omit(x)]
  
  test<- row.match(present.dat, absent.dat)
  if(length(which(is.na(test) == FALSE))) stop ('Multi-value match')
  
  #Does a scenario have that BBCH?
  present.dat$BBCH_Ref<- "Yes"
  absent.dat$BBCH_Ref<- "No"
  
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
    BBCH_Ref=as.factor(BBCH_Ref)
  )
  
  train.dat<- createDataPartition(machine.dat$BBCH_Ref, p=0.7, list=FALSE) #Select random rows to be the training data
  test.dat<- machine.dat[!train.dat] #Select other rows to be the test data
  train.dat<- machine.dat[train.dat]
  train.class<- train.dat$BBCH_Ref; train.dat<- train.dat[, -"BBCH_Ref"]
  test.class<- test.dat$BBCH_Ref; test.dat<- test.dat[, -"BBCH_Ref"]
  
  #RFE
  rf.control<- rfeControl(functions=rfFuncs, method="cv", number=10)
  rf.rfe <- rfe(train.dat, train.class,
                sizes = c(3:13),
                rfeControl = rf.control)
  
  assign(paste0("rfe_", bbch), rf.rfe)
  
  print(i)
}

#Compare differences in bbch - collect output

var.dat.out<- NULL
cm.dat.out<- NULL
imp.dat.out<- NULL
res.dat.out<- NULL

bbch.plot.list<- paste0(bbch.list, "_rfe")

pdf("220119-BBCH_Separated_RFE.pdf", onefile=TRUE)
p<- par(mfrow=c(4,3))
for(i in c(1:length(bbch.plot.list))){
  
  plot.nam<- bbch.plot.list[i]
  print(plot(get(plot.nam), type = c("o", "g"), main=paste("BBCH", gsub("_rfe", "", plot.nam))))
  
  var.dat.fill<- as.data.table(get(plot.nam)$variables)
  var.dat.fill$bbch<- plot.nam
  var.dat.out<- as.data.table(rbind(var.dat.out, var.dat.fill))

  cm.dat.fill<- as.data.table(get(plot.nam)$fit$confusion)
  cm.dat.fill$bbch<- plot.nam
  cm.dat.out<- as.data.table(rbind(cm.dat.out, cm.dat.fill))

  imp.dat.fill<- as.data.table(get(plot.nam)$optVariables)
  imp.dat.fill$bbch<- plot.nam
  imp.dat.out<- as.data.table(rbind(imp.dat.out, imp.dat.fill))

  res.dat.fill<- as.data.table(get(plot.nam)$results)
  res.dat.fill$bbch<- plot.nam
  res.dat.out<- as.data.table(rbind(res.dat.out, res.dat.fill))
  
}
par(p)
graphics.off()

setwd("C:/Users/markos/Documents/DataDownloads/BBCH_Sep_Files")
save.image(file="220119-Sep_BBCH_Analysis_RFE.RData")

setwd("C:/Users/markos/Documents/Coding_Help/Bayer_Data_Gaps")
fwrite(var.dat.out, file="220119-Full_Variables_Sep_BBCH_RFE.csv")
fwrite(cm.dat.out, file="220119-Full_CM_Sep_BBCH_RFE.csv")
fwrite(imp.dat.out, file="220119-Full_Important_Variables_Sep_BBCH_RFE.csv")
fwrite(res.dat.out, file="220119-Full_Results_Sep_BBCH_RFE.csv")

