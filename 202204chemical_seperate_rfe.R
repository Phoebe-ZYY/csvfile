############################################################################

### Chemical Separate Models RFE
###

############################################################################

rm(list = ls())

library(data.table)
library(prodlim)
library(caret)
library(randomForest)

Sys.setenv(LANG = "en")

#Test datasets

exist.dat <- fread("C:/me/data gaps/country_input_full.csv")
exist.dat$climate_main <- substr(exist.dat$layer_climate, 1, 1)
exist.dat$climate_sub <- substr(exist.dat$layer_climate, 2, 2)

exist.dat<- exist.dat[, -c("V1", "summass", "sumarea", "mindose", "name_chemclass", "wd", "name_country")]

exist.dat <- exist.dat[is.na(HDI) == FALSE & is.na(landevap) == FALSE]

chem.list<- unique(exist.dat$casrn_chemical)

#nrow(test)/nrow(bbch.dat) #work with this 84% for now
#bbch.dat<- test

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
#length(unique(bbch.dat$name_ai.x)) == length(unique(bbch.dat$casrn_chemical)) #Only use 1 of these columns

#bbch.dat<- bbch.dat[, -c("V1", "name_ai.x")]

#Check on different input data present
#table(bbch.dat$bbch_specific)

#Prep set of general scenarios that match chemical-specific scenarios
general.scenario<- exist.dat[, -c("casrn_chemical")]
general.scenario<- general.scenario[!duplicated(general.scenario)] #Lots of replicates in scenarios, remove

write.csv(exist.dat, "chemical imputation.csv")



#Now do rfe for each split and see what happens

#for(i in c(1:length(chem.list))){
for(i in c(154:300)){  
  chem<- chem.list[i]
  present.dat<- exist.dat[casrn_chemical %in% chem] #ID scenarios that use that chem
  present.dat<- present.dat[, -c("casrn_chemical")]
  present.dat<- present.dat[!duplicated(present.dat)]
  x<- row.match(present.dat, general.scenario) #ID scenarios that do not use that chem
  absent.dat<- general.scenario[!na.omit(x)]
  
  #Do any of the present rows match absent rows? If yes, there's an error
  test<- row.match(present.dat, absent.dat)
  if(length(which(is.na(test) == FALSE)) != 0) stop ('Multi-value match')
  
  #Does a scenario have that chemical?
  present.dat$ChemPresent<- "Yes"
  absent.dat$ChemPresent<- "No"
  
  machine.dat<- as.data.table(rbind(present.dat, absent.dat))
  
  machine.dat<- transform(
    machine.dat,
    #name_country=as.factor(as.character(name_country)),
    #name_crop_earthstat=as.factor(as.character(name_crop_earthstat)),
    name_cropgroup=as.factor(as.character(name_cropgroup)),
    name_region = as.factor(as.character(name_region)),
    #casrn_chemical=as.factor(as.character(casrn_chemical)),
    #class1=as.factor(as.character(class1)),
    #name_indication=as.factor(as.character(name_indication)),
    layer_climate=as.factor(as.character(layer_climate)),
    climate_sub = as.factor(as.character(climate_sub)),
    climate_main = as.factor(as.character(climate_main)),
    HDI=as.numeric(as.character(HDI)),
    #GDP_2015=as.numeric(as.character(GDP_2015)),
    landevap=as.numeric(as.character(landevap)),
    #gpw_2015=as.numeric(as.character(gpw_2015)),
    #tmin_all_mean=as.numeric(as.character(tmin_all_mean)),
    #tmax_all_mean=as.numeric(as.character(tmax_all_mean)),
    #precip_all_mean=as.numeric(as.character(precip_all_mean)),
    #ocs=as.numeric(as.character(ocs)),
    #soc_15_proj=as.numeric(as.character(soc_15_proj)),
    #soc_5_proj=as.numeric(as.character(soc_5_proj)),
    ChemPresent=as.factor(ChemPresent)
  )
  
  #machine.dat <- machine.dat %>%
    #select(ChemPresent, name_crop_earthstat, name_region, name_cropgroup,HDI, landevap, layer_climate,climate_main, climate_sub)
  
  ## make sure first column is always the label and class labels are always 1, 2, ...
  #machine.dat[, 1] = as.factor(machine.dat[, 1])
  #machine.dat[, 1] = as.numeric(unlist(machine.dat[, 1]))
  
  train.dat<- createDataPartition(machine.dat$ChemPresent, p=0.7, list=FALSE) #Select random rows to be the training data
  test.dat<- machine.dat[-train.dat] #Select other rows to be the test data
  train.dat<- machine.dat[train.dat]
  train.class<- train.dat$ChemPresent; train.dat<- train.dat[, -"ChemPresent"]
  test.class<- test.dat$ChemPresent; test.dat<- test.dat[, -"ChemPresent"]
  
  #RFE
  rf.control<- rfeControl(functions=rfFuncs, method="cv", number=10)
  rf.rfe <- rfe(train.dat, train.class,
                sizes = c(2:9),
                rfeControl = rf.control)
  
  assign(paste0("rfe_", chem), rf.rfe)
  
  print(i)
}

#Compare differences in bbch - collect output

var.dat.out<- NULL
cm.dat.out<- NULL
imp.dat.out<- NULL
res.dat.out<- NULL

chem.plot.list<- paste0(chem.list, "_rfe")
chem.plot.list <- paste0("rfe_", chem.list)

pdf("220416-chem_Separated_RFE300.pdf", onefile=TRUE)
p<- par(mfrow=c(4,3))
#for(i in c(1:length(chem.plot.list))){
for(i in c(154:300)){  
  plot.nam<- chem.plot.list[i]
  print(plot(get(plot.nam), type = c("o", "g"), main=paste("chem", gsub("_rfe", "", plot.nam))))
  
  var.dat.fill<- as.data.table(get(plot.nam)$variables)
  var.dat.fill$chem<- plot.nam
  var.dat.out<- as.data.table(rbind(var.dat.out, var.dat.fill))
  
  cm.dat.fill<- as.data.table(get(plot.nam)$fit$confusion)
  cm.dat.fill$chem<- plot.nam
  cm.dat.out<- as.data.table(rbind(cm.dat.out, cm.dat.fill))
  
  imp.dat.fill<- as.data.table(get(plot.nam)$optVariables)
  imp.dat.fill$chem<- plot.nam
  imp.dat.out<- as.data.table(rbind(imp.dat.out, imp.dat.fill))
  
  res.dat.fill<- as.data.table(get(plot.nam)$results)
  res.dat.fill$chem<- plot.nam
  res.dat.out<- as.data.table(rbind(res.dat.out, res.dat.fill))
  
}
par(p)
graphics.off()

setwd("C:/Users/markos/Documents/DataDownloads/BBCH_Sep_Files")
save.image(file="220119-Sep_BBCH_Analysis_RFE.RData")

setwd("C:/Users/markos/Documents/Coding_Help/Bayer_Data_Gaps")
fwrite(var.dat.out, file="220416-Full_Variables_Sep_BBCH_RFE300.csv")
fwrite(cm.dat.out, file="220416-Full_CM_Sep_BBCH_RFE300.csv")
fwrite(imp.dat.out, file="220416-Full_Important_Variables_Sep_BBCH_RFE300.csv")
fwrite(res.dat.out, file="220416-Full_Results_Sep_BBCH_RFE300.csv")

app.dat <- fread("C:/me/data gaps/application_input_full.csv")
app.dat<- app.dat[, -c("V1", "name_ai", "casrn_chemical")]
app.dat <- app.dat[is.na(HDI) == FALSE & is.na(landevap) == FALSE]
app.general.scenario<- app.dat[, -c("name_refapp")]
app.general.scenario<- app.general.scenario[!duplicated(app.general.scenario)]
app.present.dat<- app.dat[name_refapp %in% "Tunnel sprayer"] #ID scenarios that use that chem
app.present.dat<- app.present.dat[, -c("name_refapp")]
app.present.dat<- app.present.dat[!duplicated(app.present.dat)]
x<- row.match(app.present.dat, app.general.scenario) #ID scenarios that do not use that chem
app.absent.dat<- app.general.scenario[!na.omit(x)]


#Does a scenario have that chemical?
app.present.dat$appPresent<- "Yes"
app.absent.dat$appPresent<- "No"

app.machine.dat<- as.data.table(rbind(app.present.dat, app.absent.dat))
write.csv(app.machine.dat, "app_tunnelsprayer.csv")
unique(app.dat$name_refapp)
