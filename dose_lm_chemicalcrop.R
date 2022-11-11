############################################################################

### dose linear regression  per chemical-crop
### Yuyue Zhang

############################################################################

summary_wdhaother
annualtemp
library(dplyr)
library(data.table)

#per chemical-crop-country
lm.dat.out1<-data.frame(name_crop=c(), CASRN=c(), countcountry = c(), p_values = c(), rsquared = c(), adj_rsquared = c())
lm.dat.out<- matrix(data="X", nrow=nrow(chemcrop), ncol= 6, dimnames=list(NULL, c("name_crop", "CASRN", "countcountry", "p-values", "rsquared", "adj_rsquared")))
chemcrop <- summary_wdhaother %>%
  group_by(name_crop_earthstat, casrn_chemical) %>%
  summarize(count = n(), countcountry = n_distinct(name_country)) %>%
  filter(countcountry > 2) %>%
  select(name_crop_earthstat, casrn_chemical, countcountry) %>%
  mutate(ID = row_number())
#chemcrop[3,1]
i = 1

for(i in c(1:nrow(chemcrop))){
  crop <- chemcrop[i,1]
  chem <- chemcrop[i,2]
  countcountry <- chemcrop[i,3]
  
  lm.dat <- summary_wdhaother %>%
    #filter(wd < 2) %>%
    filter(name_crop_earthstat == crop) %>%
    filter(casrn_chemical == chem) %>%
    left_join(annualtemp, by = c("name_country" = "name_country"))# %>%
    #filter(!is.na(`annual precipitation mm`))
  
  fit <- lm(wd ~ `annual temp`  , data = lm.dat)
  
  val.out<- c(crop,
              chem, #Total Nos in dataset
              countcountry, #Total yes in dataset
              p_values = summary(fit)$coefficient[2,4],
              rsquared = summary(fit)$r.squared,
              adj_rsquared = summary(fit)$adj.r.squared
  )
  
  lm.dat.out1 <- rbind(lm.dat.out1, val.out)
  print(i)

}
summary(fit)$coefficient[2,4]
crop
write.csv(lm.dat.out1, "20220901doselinearregression_chemcrop.csv")

#per indication crop country
lm.dat.out3<-data.frame(name_crop=c(), indication=c(), countcountry = c(), p_values = c(), rsquared = c(), adj_rsquared = c())
#lm.dat.out<- matrix(data="X", nrow=nrow(chemcrop), ncol= 6, dimnames=list(NULL, c("name_crop", "CASRN", "countcountry", "p-values", "rsquared", "adj_rsquared")))
chemcrop <- summary_wdhaotheraggrindi %>%
  group_by(name_crop_earthstat, name_indication) %>%
  summarize(count = n(), countcountry = n_distinct(name_country)) %>%
  filter(countcountry > 5) %>%
  select(name_crop_earthstat, name_indication, countcountry) %>%
  mutate(ID = row_number())
#chemcrop[3,1]
i = 1

for(i in c(1:nrow(chemcrop))){
  crop <- chemcrop[i,1]
  indi <- chemcrop[i,2]
  countcountry <- chemcrop[i,3]
  
  lm.dat <- summary_wdhaotheraggrindi %>%
    #filter(wd < 2) %>%
    filter(name_crop_earthstat == crop) %>%
    filter(name_indication == indi) %>%
    left_join(annualtemp, by = c("name_country" = "name_country"))# %>%
  #filter(!is.na(`annual precipitation mm`))

  fit <- lm(wd ~ `annual temp`  , data = lm.dat)
  
  val.out<- c(crop,
              indi, 
              countcountry, 
              p_values = summary(fit)$coefficient[2,4],
              rsquared = summary(fit)$r.squared,
              adj_rsquared = summary(fit)$adj.r.squared
  )
  
  lm.dat.out3 <- rbind(lm.dat.out2, val.out)
  print(i)
  
}
write.csv(lm.dat.out3, "20220607doselinearregression_indicropannualtemp.csv")

#per chemicalclass crop country
lm.dat.out4<-data.frame(name_crop=c(), chemclass=c(), countcountry = c(), p_values = c(), rsquared = c(), adj_rsquared = c())
#lm.dat.out<- matrix(data="X", nrow=nrow(chemcrop), ncol= 6, dimnames=list(NULL, c("name_crop", "CASRN", "countcountry", "p-values", "rsquared", "adj_rsquared")))
chemcrop <- summary_wdhaother %>%
  group_by(name_crop_earthstat, name_chemclass) %>%
  summarize(count = n(), countcountry = n_distinct(name_country)) %>%
  filter(countcountry > 5) %>%
  select(name_crop_earthstat, name_chemclass, countcountry) %>%
  mutate(ID = row_number())
#chemcrop[3,1]
i = 1

for(i in c(1:nrow(chemcrop))){
  crop <- chemcrop[i,1]
  chemclass <- chemcrop[i,2]
  countcountry <- chemcrop[i,3]
  
  lm.dat <- summary_wdhaother %>%
    #filter(wd < 2) %>%
    filter(name_crop_earthstat == crop) %>%
    filter(name_chemclass == chemclass) %>%
    left_join(annualtemp, by = c("name_country" = "name_country"))# %>%
  #filter(!is.na(`annual precipitation mm`))
  
  fit <- lm(wd ~ `annual temp`  , data = lm.dat)
  
  val.out<- c(crop,
              chemclass, 
              countcountry, 
              p_values = summary(fit)$coefficient[2,4],
              rsquared = summary(fit)$r.squared,
              adj_rsquared = summary(fit)$adj.r.squared
  )
  
  lm.dat.out4 <- rbind(lm.dat.out4, val.out)
  print(i)
  
}
write.csv(lm.dat.out4, "20220609doselinearregression_chemclasscropannualtemp.csv")

#per crop group chemical country
lm.dat.out5<-data.frame(chem = c(),name_cropgroup=c(), count = c(),countcountry = c(), p_values = c(), rsquared = c(), adj_rsquared = c())
#lm.dat.out<- matrix(data="X", nrow=nrow(chemcrop), ncol= 6, dimnames=list(NULL, c("name_crop", "CASRN", "countcountry", "p-values", "rsquared", "adj_rsquared")))
chemcropgroup <- summary_wdhaother %>%
  group_by(casrn_chemical, name_cropgroup,) %>%
  summarize(count = n(), countcountry = n_distinct(name_country)) %>%
  filter(count > 5) %>%
  select(casrn_chemical, name_cropgroup, count, countcountry) %>%
  mutate(ID = row_number())

i = 1

for(i in c(1:nrow(chemcropgroup))){
  chem <- chemcropgroup[i,1]
  cropgroup <- chemcropgroup[i,2]
  count <- chemcropgroup[i,3]
  countcountry <- chemcropgroup[i,4]
  
  lm.dat <- summary_wdhaother %>%
    #filter(wd < 2) %>%
    filter(casrn_chemical == chem) %>%
    filter(name_cropgroup == cropgroup) %>%
    left_join(annualtemp, by = c("name_country" = "name_country"))# %>%
  #filter(!is.na(`annual precipitation mm`))
  
  fit <- lm(wd ~ `annual temp`  , data = lm.dat)
  
  val.out<- c(chem,
              cropgroup, 
              count,
              countcountry, 
              p_values = summary(fit)$coefficient[2,4],
              rsquared = summary(fit)$r.squared,
              adj_rsquared = summary(fit)$adj.r.squared
  )
  
  lm.dat.out5 <- rbind(lm.dat.out5, val.out)
  print(i)
  
}
write.csv(lm.dat.out5, "20220620doselinearregression_chemcropgroupannualtemp.csv")

#per crop-chemical-country 95 percentile
summary_95percentile <- input_wdhaother %>%
  group_by(name_country, name_ai, name_chemclass, name_crop_earthstat, name_cropgroup, casrn_chemical, 
           name_indication) %>%
  summarize(quantile = quantile(dose, probs = 0.95))

lm.dat.out7<-data.frame(name_crop=c(), CASRN=c(), countcountry = c(), p_values = c(), rsquared = c(), adj_rsquared = c())
#lm.dat.out<- matrix(data="X", nrow=nrow(chemcrop), ncol= 6, dimnames=list(NULL, c("name_crop", "CASRN", "countcountry", "p-values", "rsquared", "adj_rsquared")))
chemcrop <- summary_95percentile %>%
  group_by(name_crop_earthstat, casrn_chemical) %>%
  summarize(count = n(), countcountry = n_distinct(name_country)) %>%
  filter(countcountry > 2) %>%
  select(name_crop_earthstat, casrn_chemical, countcountry) %>%
  mutate(ID = row_number())
#chemcrop[3,1]
i = 1

for(i in c(698:nrow(chemcrop))){
  crop <- chemcrop[i,1]
  chem <- chemcrop[i,2]
  countcountry <- chemcrop[i,3]
  
  lm.dat <- summary_95percentile %>%
    #filter(wd < 2) %>%
    filter(name_crop_earthstat == crop) %>%
    filter(casrn_chemical == chem) %>%
    #filter(name_indication == chem) %>%
    left_join(annualtemp, by = c("name_country" = "name_country"))# %>%
  #filter(!is.na(`annual precipitation mm`))
  
  fit <- lm(quantile ~ `minimum temp`  , data = lm.dat)
  
  val.out<- c(crop,
              chem, #Total Nos in dataset
              countcountry, #Total yes in dataset
              p_values = summary(fit)$coefficient[2,4],
              rsquared = summary(fit)$r.squared,
              adj_rsquared = summary(fit)$adj.r.squared
  )
  
  lm.dat.out7 <- rbind(lm.dat.out7, val.out)
  print(i)
  
}
summary(fit)$coefficient[2,4]
crop
write.csv(lm.dat.out7, "20220627doselinearregression95percentile_chemcrop.csv")



#per crop group chemclass, per group a dose, with latitude/temp/precipitation/hdi/gdp
lm.dat.out9<-data.frame(cropgroup=c(), chemclass=c(), countcountry = c(), count = c(), quadraticp_values = c(),
                        quadraticrsquared = c(),
                        quadraticadj_rsquared = c(),
                        annualtempp_values = c(),
                        annualtemprsquared = c(),
                        annualtempadj_rsquared = c(),
                        mintempp_values =c(),
                        mintemprsquared = c(),
                        mintempadj_rsquared = c(),
                        precipp_values = c(),
                        preciprsquared = c(),
                        precipadj_rsquared = c(),
                        gdpp_values = c(),
                        gdprsquared = c(),
                        gdpadj_rsquared = c(),
                        gdppercapitap_values = c(),
                        gdppercapitarsquared = c(),
                        gdppercapitaadj_rsquared = c(),
                        hdip_values = c(),
                        hdirsquared = c(),
                        hdiadj_rsquared = c())

chemcrop <- re_cropchemical %>%
  group_by(name_cropgroup, name_chemclass1) %>%
  summarize(count = n(), countcountry = n_distinct(name_country)) %>%
  filter(countcountry > 2) %>%
  select(name_cropgroup, name_chemclass1, countcountry, count) %>%
  mutate(ID = row_number())

i <- 50

for(i in c(50:nrow(chemcrop))){
  cropgroup <- chemcrop[i,1]
  chemclass <- chemcrop[i,2]
  countcountry <- chemcrop[i,3]
  count <- chemcrop[i,4]
  
  summary_95percentile <- re_cropchemical %>%
    filter(name_cropgroup == cropgroup$name_cropgroup) %>%
    filter(name_chemclass1 == chemclass$name_chemclass1) #%>%
    #group_by(name_country) %>%
    #filter(wd <= quantile(wd, 0.95) & wd >= quantile(wd, 0.05))
  
  summary_95percentile$latitude2 <- summary_95percentile$latitude ^2
  summary_95percentile$logwd <- log10(summary_95percentile$wd)
  
  quadraticModel <- lm(logwd ~ latitude + latitude2, data=summary_95percentile)
  annualtempModel <- lm(logwd ~ `annual temp`, data=summary_95percentile)
  mintempModel <- lm(logwd ~ `minimum temp`, data=summary_95percentile)
  precipModel <- lm(logwd ~ `annual precipitation mm`, data=summary_95percentile)
  gdpModel <- lm(logwd ~ log(`gdp`), data=summary_95percentile)
  gdppercapitaModel <- lm(logwd ~ gdppercapita, data=summary_95percentile)
  hdiModel <- lm(logwd ~ hdi, data=summary_95percentile)
  
  val.out<- c(cropgroup,
              chemclass, 
              countcountry,
              count,
              quadraticp_values = summary(quadraticModel)$coefficient[2,4],
              quadraticrsquared = summary(quadraticModel)$r.squared,
              quadraticadj_rsquared = summary(quadraticModel)$adj.r.squared,
              annualtempp_values = summary(annualtempModel)$coefficient[2,4],
              annualtemprsquared = summary(annualtempModel)$r.squared,
              annualtempadj_rsquared = summary(annualtempModel)$adj.r.squared,
              mintempp_values = summary(mintempModel)$coefficient[2,4],
              mintemprsquared = summary(mintempModel)$r.squared,
              mintempadj_rsquared = summary(mintempModel)$adj.r.squared,
              precipp_values = summary(precipModel)$coefficient[2,4],
              preciprsquared = summary(precipModel)$r.squared,
              precipadj_rsquared = summary(precipModel)$adj.r.squared,
              gdpp_values = summary(gdpModel)$coefficient[2,4],
              gdprsquared = summary(gdpModel)$r.squared,
              gdpadj_rsquared = summary(gdpModel)$adj.r.squared,
              gdppercapitap_values = summary(gdppercapitaModel)$coefficient[2,4],
              gdppercapitarsquared = summary(gdppercapitaModel)$r.squared,
              gdppercapitaadj_rsquared = summary(gdppercapitaModel)$adj.r.squared,
              hdip_values = summary(hdiModel)$coefficient[2,4],
              hdirsquared = summary(hdiModel)$r.squared,
              hdiadj_rsquared = summary(hdiModel)$adj.r.squared
  )
  
  
  lm.dat.out9 <- rbind(lm.dat.out9, val.out)
  print(i)
}

write.csv(lm.dat.out8, "20220902cropgroupchemclass.csv")

#per crop group chemclass, 95%, with latitude/temp/precipitation/hdi/gdp
lm.dat.out8<-data.frame(cropgroup=c(), chemclass=c(), countcountry = c(), count = c(), quadraticp_values = c(),
                        quadraticrsquared = c(),
                        quadraticadj_rsquared = c(),
                        annualtempp_values = c(),
                        annualtemprsquared = c(),
                        annualtempadj_rsquared = c(),
                        mintempp_values =c(),
                        mintemprsquared = c(),
                        mintempadj_rsquared = c(),
                        precipp_values = c(),
                        preciprsquared = c(),
                        precipadj_rsquared = c(),
                        gdpp_values = c(),
                        gdprsquared = c(),
                        gdpadj_rsquared = c(),
                        gdppercapitap_values = c(),
                        gdppercapitarsquared = c(),
                        gdppercapitaadj_rsquared = c(),
                        hdip_values = c(),
                        hdirsquared = c(),
                        hdiadj_rsquared = c())

chemcrop <- re_cropchemical %>%
  group_by(name_cropgroup, name_chemclass) %>%
  summarize(count = n(), countcountry = n_distinct(name_country)) %>%
  filter(countcountry > 2) %>%
  select(name_cropgroup, name_chemclass, countcountry, count) %>%
  mutate(ID = row_number())

i <- 1

for(i in c(465:nrow(chemcrop))){
  cropgroup <- chemcrop[i,1]
  chemclass <- chemcrop[i,2]
  countcountry <- chemcrop[i,3]
  count <- chemcrop[i,4]
  
  summary_95percentile <- re_cropchemical %>%
    filter(name_cropgroup == cropgroup) %>%
    filter(name_chemclass == chemclass) #%>%
  #group_by(name_country) %>%
  #filter(wd <= quantile(wd, 0.95) & wd >= quantile(wd, 0.05))
  
  summary_95percentile$latitude2 <- summary_95percentile$latitude ^2
  summary_95percentile$logwd <- log10(summary_95percentile$wd)
  
  quadraticModel <- lm(logwd ~ latitude + latitude2, data=summary_95percentile)
  annualtempModel <- lm(logwd ~ `annual temp`, data=summary_95percentile)
  mintempModel <- lm(logwd ~ log10(`minimum temp`), data=summary_95percentile)
  precipModel <- lm(logwd ~ `annual precipitation mm`, data=summary_95percentile)
  gdpModel <- lm(logwd ~ log(`gdp`), data=summary_95percentile)
  gdppercapitaModel <- lm(logwd ~ gdppercapita, data=summary_95percentile)
  hdiModel <- lm(logwd ~ hdi, data=summary_95percentile)
  
  val.out<- c(cropgroup,
              chemclass, 
              countcountry,
              count,
              quadraticp_values = summary(quadraticModel)$coefficient[2,4],
              quadraticrsquared = summary(quadraticModel)$r.squared,
              quadraticadj_rsquared = summary(quadraticModel)$adj.r.squared,
              annualtempp_values = summary(annualtempModel)$coefficient[2,4],
              annualtemprsquared = summary(annualtempModel)$r.squared,
              annualtempadj_rsquared = summary(annualtempModel)$adj.r.squared,
              mintempp_values = summary(mintempModel)$coefficient[2,4],
              mintemprsquared = summary(mintempModel)$r.squared,
              mintempadj_rsquared = summary(mintempModel)$adj.r.squared,
              precipp_values = summary(precipModel)$coefficient[2,4],
              preciprsquared = summary(precipModel)$r.squared,
              precipadj_rsquared = summary(precipModel)$adj.r.squared,
              gdpp_values = summary(gdpModel)$coefficient[2,4],
              gdprsquared = summary(gdpModel)$r.squared,
              gdpadj_rsquared = summary(gdpModel)$adj.r.squared,
              gdppercapitap_values = summary(gdppercapitaModel)$coefficient[2,4],
              gdppercapitarsquared = summary(gdppercapitaModel)$r.squared,
              gdppercapitaadj_rsquared = summary(gdppercapitaModel)$adj.r.squared,
              hdip_values = summary(hdiModel)$coefficient[2,4],
              hdirsquared = summary(hdiModel)$r.squared,
              hdiadj_rsquared = summary(hdiModel)$adj.r.squared
  )
  
  
  lm.dat.out8 <- rbind(lm.dat.out8, val.out)
  print(i)
}

write.csv(lm.dat.out8, "20220902cropgroupchemclass.csv")

