#--------------------------------------------------------------
#Ben Neely
#04/25/2025
#Investigate what drives angler decision to harvest legal and illegal fish from creel data
#This incorporates creel data from 1997-2024 in Kansas waters
#--------------------------------------------------------------

## Clear R
cat("\014")  
rm(list=ls())

## Install and load packages
## Checks if package is installed, installs if not, activates for current session
if("FSA" %in% rownames(installed.packages()) == FALSE) {install.packages("FSA")}
library(FSA)

if("rio" %in% rownames(installed.packages()) == FALSE) {install.packages("rio")}
library(rio)

if("patchwork" %in% rownames(installed.packages()) == FALSE) {install.packages("patchwork")}
library(patchwork)

if("rpart" %in% rownames(installed.packages()) == FALSE) {install.packages("rpart")}
library(rpart)

if("rpart.plot" %in% rownames(installed.packages()) == FALSE) {install.packages("rpart.plot")}
library(rpart.plot)

if("tidyverse" %in% rownames(installed.packages()) == FALSE) {install.packages("tidyverse")}
library(tidyverse)

## Set seed for reproducible results
addTaskCallback(function(...) {set.seed(1114);TRUE})

## Set ggplot theme
pubtheme=theme_classic()+
  theme(panel.grid=element_blank(), 
        panel.background=element_blank(),
        plot.background=element_blank(),
        panel.border=element_rect(color="black",fill="transparent"),
        axis.title=element_text(size=22,color="black",face="bold"),
        axis.text=element_text(size=18,color="black"),
        axis.line=element_line(color="black"),
        axis.ticks=element_line(color="black"),
        legend.position="none")
options(scipen=999)

## Set working directory
setwd("C:/Users/Ben.Neely/OneDrive - State of Kansas, OITS/Desktop/Harvest patterns/")

## Read in data with import
dat=import("harvest factors/moddat_clean.csv")

## Manipulate variables to make analyses easier
dat1=dat%>%
  mutate(harv=factor(harv),
         sppgrp=factor(sppgrp),
         legal=factor(legal),
         target_spp=factor(target_spp),
         generalist=factor(generalist),
         dayofweek=factor(dayofweek),
         day_period=factor(day_period),
         bs=factor(bs),
         hrs=as.numeric(hrs),
         instate=factor(instate))

################################################################################
################################################################################
## Predict harvest based on other variables
## We set maxdepth to five meaning there can't be more than six layers
## We set minbucket to 1% of observations in the training data set
## A node can't be split unless it has at least 1% of the number of obs in the training set
###############################################################
## Global model
## Separate into training (70%) and prediction (30%) data sets
dat1$id=1:nrow(dat1)
train=dat1%>%sample_frac(0.70)
test=anti_join(dat1,train,by="id")

## Use rpart to build recursive partitioning tree
mod=rpart(harv~sppgrp+tl+target_spp+generalist+legal+doy+dayofweek+day_period+year+
            bs+hrs+instate+tot_ang+prop_male+prop_o65+dist_km+ha+long+lat+
            pop_2010,
          parms=list(split="gini"),
          data=train,
          method="class",
          maxdepth=5,
          minbucket=nrow(train)*0.01)

## Look at classification tree data
## https://www.r-bloggers.com/2022/10/understanding-leaf-node-numbers-when-using-rpart-and-rpart-rules/
## http://www.milbo.org/rpart-plot/prp.pdf
## Extra = 2 shows the number of individuals with the title and the total number of individuals in the node
## i.e., if the box says release, it shows number of fish released on the left
## Dividing these numbers gives proportion released or harvested (title of node)
rpart.plot(mod,extra=2)
rpart.plot(mod,extra=8)
rpart.rules(mod)

## Fit model to test data and look at correct classification rate
test$pred=predict(mod,test,type="class")
mean(test$harv==test$pred)
nrow(dat1)

## Variable importance and relative variable importance
imp_g=tibble(var=names(mod$variable.importance),vi=mod$variable.importance)%>%
  mutate(prop_vi=vi/sum(vi),
         rel_prop_vi=vi/max(vi))

###############################################################
###############################################################
###############################################################
## Combine and examine variable importance
## List of all possible variables
varlist=c("tl","target_spp","generalist","legal","doy","dayofweek",
          "day_period","year","bs","hrs","instate","tot_ang","prop_male",
          "prop_o65","dist_km","ha","long","lat","pop_2010")

## Add all variables to all species and make value 0 if not included in model
varimp_g=imp_g%>%
  complete(var=varlist,
           fill=list(vi=0,
                     prop_vi=0,
                     rel_prop_vi=0))

## Create plotting data
varimp_plotdat_g=varimp_g%>%
  mutate(name=fct_reorder(var,rel_prop_vi),
         name=recode(name,
                     "sppgrp"="Species",
                     "tl"="Length",
                     "target_spp"="Target species",
                     "generalist"="Generalist",
                     "legal"="Legal",
                     "doy"="Day of year",
                     "dayofweek"="Day of week",
                     "day_period"="Time of day",
                     "year"="Year",
                     "bs"="Boat or shore",
                     "hrs"="Hours angling",
                     "instate"="In-state",
                     "tot_ang"="Number in party",
                     "prop_male"="Proportion male",
                     "prop_o65"="Proportion over 65",
                     "dist_km"="Distance traveled",
                     "ha"="Surface area",
                     "long"="Longitude",
                     "lat"="Latitude",
                     "pop_2010"="County population"),
         spp=case_when(name=="Species"~"spp", TRUE~"no"))

## Create global variable importance plot
a=ggplot(varimp_plotdat_g)+
  geom_bar(aes(x=prop_vi,y=name),stat="identity",fill="gray30")+
  scale_x_continuous(breaks=seq(0,0.5,0.1))+
  coord_cartesian(xlim=c(-0.005,0.52),
                  ylim=c(0.4,20.6),
                  expand=F)+
  labs(x="Variable importance",
       y="")+
  pubtheme

############################################################################################
############################################################################################
## Set up sppgrp specific data sets
bcf=filter(dat1,sppgrp=="bcf")
ccarp=filter(dat1,sppgrp=="ccarp")
ccf=filter(dat1,sppgrp=="ccf")
crappie=filter(dat1,sppgrp=="crappie")
fhc=filter(dat1,sppgrp=="fhc")
fwd=filter(dat1,sppgrp=="fwd")
lmb=filter(dat1,sppgrp=="lmb")
percid=filter(dat1,sppgrp=="percids")
smb=filter(dat1,sppgrp=="smb")
sunfish=filter(dat1,sppgrp=="sunfish")
whb=filter(dat1,sppgrp=="whb")
wiper=filter(dat1,sppgrp=="wiper")

################################################################################
################################################################################
## Predict harvest based on other variables
## We set maxdepth to five meaning there can't be more than six layers
## We set minbucket to 1% of the total observations in the training data set
###############################################################
## Blue Catfish
## Separate into training (70%) and prediction (30%) data sets
bcf$id=1:nrow(bcf)
bcftrain=bcf%>%sample_frac(0.70)
bcftest=anti_join(bcf,bcftrain,by="id")

## Use rpart to build recursive partitioning tree
bcftree=rpart(harv~tl+legal+target_spp+generalist+doy+dayofweek+day_period+year+
                bs+hrs+instate+tot_ang+prop_male+prop_o65+dist_km+ha+long+lat+
                pop_2010,
              parms=list(split="gini"),
              data=bcftrain,
              method="class",
              maxdepth=5,
              minbucket=nrow(bcftrain)*0.01)

## Look at classification tree data
## https://www.r-bloggers.com/2022/10/understanding-leaf-node-numbers-when-using-rpart-and-rpart-rules/
## http://www.milbo.org/rpart-plot/prp.pdf
## Extra = 2 shows the number of individuals with the title and the total number of individuals in the node
## i.e., if the box says release, it shows number of fish released on the left
## Dividing these numbers gives proportion released or harvested (title of node)
## This is also showning the extra = 8 tree
rpart.plot(bcftree,extra=2)
rpart.plot(bcftree,extra=8)
rpart.rules(bcftree)

## Fit model to test data and look at correct classification rate
bcftest$pred=predict(bcftree,bcftest,type="class")
mean(bcftest$harv==bcftest$pred)
nrow(bcftrain)

## Variable importance and relative variable importance
imp_bcf=tibble(var=names(bcftree$variable.importance),vi=bcftree$variable.importance)%>%
  mutate(spp="Blue Catfish",.before=1,
         prop_vi=vi/sum(vi),
         rel_prop_vi=vi/max(vi))%>%
  select(spp,var,vi,prop_vi,rel_prop_vi)

###############################################################
## Channel Catfish
## Separate into training (70%) and prediction (30%) data sets
ccf$id=1:nrow(ccf)
ccftrain=ccf%>%sample_frac(0.70)
ccftest=anti_join(ccf,ccftrain,by="id")

## Use rpart to build recursive partitioning tree
ccftree=rpart(harv~tl+legal+target_spp+generalist+doy+dayofweek+day_period+year+
                bs+hrs+instate+tot_ang+prop_male+prop_o65+dist_km+ha+long+lat+
                pop_2010,
              parms=list(split="gini"),
              data=ccftrain,
              method="class",
              maxdepth=5,
              minbucket=nrow(ccftrain)*0.01)

## Look at classification tree data
rpart.plot(ccftree,extra=2)
rpart.plot(ccftree,extra=8)
rpart.rules(ccftree)

## Fit model to test data and look at correct classification rate
ccftest$pred=predict(ccftree,ccftest,type="class")
mean(ccftest$harv==ccftest$pred)
nrow(ccftrain)

## Variable importance and relative variable importance
imp_ccf=tibble(var=names(ccftree$variable.importance),vi=ccftree$variable.importance)%>%
  mutate(spp="Channel Catfish",.before=1,
         prop_vi=vi/sum(vi),
         rel_prop_vi=vi/max(vi))%>%
  select(spp,var,vi,prop_vi,rel_prop_vi)

###############################################################
## Common Carp
## Separate into training (70%) and prediction (30%) data sets
ccarp$id=1:nrow(ccarp)
ccarptrain=ccarp%>%sample_frac(0.70)
ccarptest=anti_join(ccarp,ccarptrain,by="id")

## Use rpart to build recursive partitioning tree
ccarptree=rpart(harv~tl+legal+target_spp+generalist+doy+dayofweek+day_period+year+
                  bs+hrs+instate+tot_ang+prop_male+prop_o65+dist_km+ha+long+lat+
                  pop_2010,
                parms=list(split="gini"),
                data=ccarptrain,
                method="class",
                maxdepth=5,
                minbucket=nrow(ccarptrain)*0.01)

## Look at classification tree data
rpart.plot(ccarptree,extra=2)
rpart.plot(ccarptree,extra=8)
rpart.rules(ccarptree)

## Fit model to test data and look at correct classification rate
ccarptest$pred=predict(ccarptree,ccarptest,type="class")
mean(ccarptest$harv==ccarptest$pred)
nrow(ccarptrain)

## Variable importance and relative variable importance
## Divide all variable importance scores by the max observed
imp_ccarp=tibble(var=names(ccarptree$variable.importance),vi=ccarptree$variable.importance)%>%
  mutate(spp="Common Carp",.before=1,
         prop_vi=vi/sum(vi),
         rel_prop_vi=vi/max(vi))%>%
  select(spp,var,vi,prop_vi,rel_prop_vi)

###############################################################
## Crappies
## Separate into training (70%) and prediction (30%) data sets
crappie$id=1:nrow(crappie)
crappietrain=crappie%>%sample_frac(0.70)
crappietest=anti_join(crappie,crappietrain,by="id")

## Use rpart to build recursive partitioning tree
crappietree=rpart(harv~tl+legal+target_spp+generalist+doy+dayofweek+day_period+year+
                bs+hrs+instate+tot_ang+prop_male+prop_o65+dist_km+ha+long+lat+
                pop_2010,
              parms=list(split="gini"),
              data=crappietrain,
              method="class",
              maxdepth=5,
              minbucket=nrow(crappietrain)*0.01)

## Look at classification tree data
rpart.plot(crappietree,extra=2)
rpart.plot(crappietree,extra=8)
rpart.rules(crappietree)

## Fit model to test data and look at correct classification rate
crappietest$pred=predict(crappietree,crappietest,type="class")
mean(crappietest$harv==crappietest$pred)
nrow(crappietrain)

## Variable importance and relative variable importance
imp_crappie=tibble(var=names(crappietree$variable.importance),vi=crappietree$variable.importance)%>%
  mutate(spp="Crappies",.before=1,
         prop_vi=vi/sum(vi),
         rel_prop_vi=vi/max(vi))%>%
  select(spp,var,vi,prop_vi,rel_prop_vi)

###############################################################
## Flathead Catfish
## Separate into training (70%) and prediction (30%) data sets
fhc$id=1:nrow(fhc)
fhctrain=fhc%>%sample_frac(0.70)
fhctest=anti_join(fhc,fhctrain,by="id")

## Use rpart to build recursive partitioning tree
fhctree=rpart(harv~tl+legal+target_spp+generalist+doy+dayofweek+day_period+year+
                bs+hrs+instate+tot_ang+prop_male+prop_o65+dist_km+ha+long+lat+
                pop_2010,
              parms=list(split="gini"),
              data=fhctrain,
              method="class",
              maxdepth=5,
              minbucket=nrow(fhctrain)*0.01)

## Look at classification tree data
rpart.plot(fhctree,extra=2)
rpart.plot(fhctree,extra=8)
rpart.rules(fhctree)

## Fit model to test data and look at correct classification rate
fhctest$pred=predict(fhctree,fhctest,type="class")
mean(fhctest$harv==fhctest$pred)
nrow(fhctrain)

## Variable importance and relative variable importance
imp_fhc=tibble(var=names(fhctree$variable.importance),vi=fhctree$variable.importance)%>%
  mutate(spp="Flathead Catfish",.before=1,
         prop_vi=vi/sum(vi),
         rel_prop_vi=vi/max(vi))%>%
  select(spp,var,vi,prop_vi,rel_prop_vi)

###############################################################
## Freshwater Drum
## Separate into training (70%) and prediction (30%) data sets
fwd$id=1:nrow(fwd)
fwdtrain=fwd%>%sample_frac(0.70)
fwdtest=anti_join(fwd,fwdtrain,by="id")

## Use rpart to build recursive partitioning tree
fwdtree=rpart(harv~tl+legal+target_spp+generalist+doy+dayofweek+day_period+year+
                bs+hrs+instate+tot_ang+prop_male+prop_o65+dist_km+ha+long+lat+
                pop_2010,
              parms=list(split="gini"),
              data=fwdtrain,
              method="class",
              maxdepth=5,
              minbucket=nrow(fwdtrain)*0.01)

## Look at classification tree data
rpart.plot(fwdtree,extra=2)
rpart.plot(fwdtree,extra=8)
rpart.rules(fwdtree)

## Fit model to test data and look at correct classification rate
fwdtest$pred=predict(fwdtree,fwdtest,type="class")
mean(fwdtest$harv==fwdtest$pred)
nrow(fwdtrain)

## Variable importance and relative variable importance
imp_fwd=tibble(var=names(fwdtree$variable.importance),vi=fwdtree$variable.importance)%>%
  mutate(spp="Freshwater Drum",.before=1,
         prop_vi=vi/sum(vi),
         rel_prop_vi=vi/max(vi))%>%
  select(spp,var,vi,prop_vi,rel_prop_vi)

###############################################################
## Largemouth Bass
## Separate into training (70%) and prediction (30%) data sets
lmb$id=1:nrow(lmb)
lmbtrain=lmb%>%sample_frac(0.70)
lmbtest=anti_join(lmb,lmbtrain,by="id")

## Use rpart to build recursive partitioning tree
lmbtree=rpart(harv~tl+legal+target_spp+generalist+doy+dayofweek+day_period+year+
                bs+hrs+instate+tot_ang+prop_male+prop_o65+dist_km+ha+long+lat+
                pop_2010,
              parms=list(split="gini"),
              data=lmbtrain,
              method="class",
              maxdepth=5,
              minbucket=nrow(lmbtrain)*0.01)

## Look at classification tree data
rpart.plot(lmbtree,extra=2)
rpart.plot(lmbtree,extra=8)
rpart.rules(lmbtree)

## Fit model to test data and look at correct classification rate
lmbtest$pred=predict(lmbtree,lmbtest,type="class")
mean(lmbtest$harv==lmbtest$pred)
nrow(lmbtrain)

## No variables to gauge importance

###############################################################
## Palmetto Bass
## Separate into training (70%) and prediction (30%) data sets
wiper$id=1:nrow(wiper)
wipertrain=wiper%>%sample_frac(0.70)
wipertest=anti_join(wiper,wipertrain,by="id")

## Use rpart to build recursive partitioning tree
wipertree=rpart(harv~tl+legal+target_spp+generalist+doy+dayofweek+day_period+year+
                  bs+hrs+instate+tot_ang+prop_male+prop_o65+dist_km+ha+long+lat+
                  pop_2010,
                parms=list(split="gini"),
                data=wipertrain,
                method="class",
                maxdepth=5,
                minbucket=nrow(wipertrain)*0.01)

## Look at classification tree data
rpart.plot(wipertree,extra=2)
rpart.plot(wipertree,extra=8)
rpart.rules(wipertree)

## Fit model to test data and look at correct classification rate
wipertest$pred=predict(wipertree,wipertest,type="class")
mean(wipertest$harv==wipertest$pred)
nrow(wipertrain)

## Variable importance and relative variable importance
imp_wiper=tibble(var=names(wipertree$variable.importance),vi=wipertree$variable.importance)%>%
  mutate(spp="Palmetto Bass",.before=1,
         prop_vi=vi/sum(vi),
         rel_prop_vi=vi/max(vi))%>%
  select(spp,var,vi,prop_vi,rel_prop_vi)


###############################################################
## Percids
## Separate into training (70%) and prediction (30%) data sets
percid$id=1:nrow(percid)
percidtrain=percid%>%sample_frac(0.70)
percidtest=anti_join(percid,percidtrain,by="id")

## Use rpart to build recursive partitioning tree
percidtree=rpart(harv~tl+legal+target_spp+generalist+doy+dayofweek+day_period+year+
                  bs+hrs+instate+tot_ang+prop_male+prop_o65+dist_km+ha+long+lat+
                  pop_2010,
                parms=list(split="gini"),
                data=percidtrain,
                method="class",
                maxdepth=5,
                minbucket=nrow(percidtrain)*0.01)

## Look at classification tree data
rpart.plot(percidtree,extra=2)
rpart.plot(percidtree,extra=8)
rpart.rules(percidtree)

## Fit model to test data and look at correct classification rate
percidtest$pred=predict(percidtree,percidtest,type="class")
mean(percidtest$harv==percidtest$pred)
nrow(percidtrain)

## Variable importance and relative variable importance
imp_percid=tibble(var=names(percidtree$variable.importance),vi=percidtree$variable.importance)%>%
  mutate(spp="Percids",.before=1,
         prop_vi=vi/sum(vi),
         rel_prop_vi=vi/max(vi))%>%
  select(spp,var,vi,prop_vi,rel_prop_vi)

###############################################################
## Smallmouth Bass
## Separate into training (70%) and prediction (30%) data sets
smb$id=1:nrow(smb)
smbtrain=smb%>%sample_frac(0.70)
smbtest=anti_join(smb,smbtrain,by="id")

## Use rpart to build recursive partitioning tree
smbtree=rpart(harv~tl+legal+target_spp+generalist+doy+dayofweek+day_period+year+
                bs+hrs+instate+tot_ang+prop_male+prop_o65+dist_km+ha+long+lat+
                pop_2010,
              parms=list(split="gini"),
              data=smbtrain,
              method="class",
              maxdepth=5,
              minbucket=nrow(smbtrain)*0.01)

## Look at classification tree data
rpart.plot(smbtree,extra=2)
rpart.plot(smbtree,extra=8)
rpart.rules(smbtree)

## Fit model to test data and look at correct classification rate
smbtest$pred=predict(smbtree,smbtest,type="class")
mean(smbtest$harv==smbtest$pred)
nrow(smbtrain)

## No variables to gauge importance

###############################################################
## Sunfishes
## Separate into training (70%) and prediction (30%) data sets
sunfish$id=1:nrow(sunfish)
sunfishtrain=sunfish%>%sample_frac(0.70)
sunfishtest=anti_join(sunfish,sunfishtrain,by="id")

## Use rpart to build recursive partitioning tree
sunfishtree=rpart(harv~tl+legal+target_spp+generalist+doy+dayofweek+day_period+year+
                    bs+hrs+instate+tot_ang+prop_male+prop_o65+dist_km+ha+long+lat+
                    pop_2010,
                  parms=list(split="gini"),
                  data=sunfishtrain,
                  method="class",
                  maxdepth=5,
                  minbucket=nrow(sunfishtrain)*0.01)

## Look at classification tree data
rpart.plot(sunfishtree,extra=2)
rpart.plot(sunfishtree,extra=8)
rpart.rules(sunfishtree)

## Fit model to test data and look at correct classification rate
sunfishtest$pred=predict(sunfishtree,sunfishtest,type="class")
mean(sunfishtest$harv==sunfishtest$pred)
nrow(sunfishtrain)

## Variable importance and relative variable importance
imp_sunfish=tibble(var=names(sunfishtree$variable.importance),vi=sunfishtree$variable.importance)%>%
  mutate(spp="Sunfishes",.before=1,
         prop_vi=vi/sum(vi),
         rel_prop_vi=vi/max(vi))%>%
  select(spp,var,vi,prop_vi,rel_prop_vi)

###############################################################
## White Bass
## Separate into training (70%) and prediction (30%) data sets
whb$id=1:nrow(whb)
whbtrain=whb%>%sample_frac(0.70)
whbtest=anti_join(whb,whbtrain,by="id")

## Use rpart to build recursive partitioning tree
whbtree=rpart(harv~tl+legal+target_spp+generalist+doy+dayofweek+day_period+year+
                bs+hrs+instate+tot_ang+prop_male+prop_o65+dist_km+ha+long+lat+
                pop_2010,
              parms=list(split="gini"),
              data=whbtrain,
              method="class",
              maxdepth=5,
              minbucket=nrow(whbtrain)*0.01)

## Look at classification tree data
rpart.plot(whbtree,extra=2)
rpart.plot(whbtree,extra=8)
rpart.rules(whbtree)

## Fit model to test data and look at correct classification rate
whbtest$pred=predict(whbtree,whbtest,type="class")
mean(whbtest$harv==whbtest$pred)
nrow(whbtrain)

## Variable importance and relative variable importance
imp_whb=tibble(var=names(whbtree$variable.importance),vi=whbtree$variable.importance)%>%
  mutate(spp="White Bass",.before=1,
         prop_vi=vi/sum(vi),
         rel_prop_vi=vi/max(vi))%>%
  select(spp,var,vi,prop_vi,rel_prop_vi)

###############################################################
###############################################################
###############################################################
## Combine and examine variable importance for each species/taxa model
## Set up a list of all species/species groups
spplist=c("Blue Catfish","Channel Catfish","Common Carp","Crappies",
          "Flathead Catfish","Freshwater Drum","Largemouth Bass",
          "Palmetto Bass","Percids","Smallmouth Bass","Sunfishes","White Bass")

## List of all possible variables
varlist=c("tl","target_spp","generalist","legal","doy","dayofweek",
          "day_period","year","bs","hrs","instate","tot_ang","prop_male",
          "prop_o65","dist_km","ha","long","lat","pop_2010")

## Add all variables to all species and make value 0 if not included in model
varimp=bind_rows(imp_bcf,imp_ccarp,imp_ccf,imp_crappie,
                 imp_fhc,imp_fwd,imp_percid,
                 imp_sunfish,imp_whb,imp_wiper)%>%
  complete(spp=spplist,
           var=varlist,
           fill=list(vi=0,
                     prop_vi=0,
                     rel_prop_vi=0))

################################################################################
## Create plotting data
## Calculating mean variable importance across the 12 models with 95% CI
varimp_plotdat=varimp%>%
  group_by(var)%>%
  summarize(mods=n(),
            mean_imp=mean(prop_vi),
            sd_imp=sd(prop_vi),
            lci95=mean_imp-(qt(0.975,df=mods-1)*sd_imp/sqrt(mods)),
            uci95=mean_imp+(qt(0.975,df=mods-1)*sd_imp/sqrt(mods)))%>%
  ungroup()%>%
  mutate(name=fct_reorder(var,mean_imp),
         name=recode(name,
                     "sppgrp"="Species",
                     "tl"="Length (10)",
                     "target_spp"="Target species (3)",
                     "generalist"="Generalist (1)",
                     "legal"="Legal (5)",
                     "doy"="Day of year (8)",
                     "dayofweek"="Day of week (0)",
                     "day_period"="Time of day (7)",
                     "year"="Year (7)",
                     "bs"="Boat or shore (3)",
                     "hrs"="Hours angling (5)",
                     "instate"="In-state (3)",
                     "tot_ang"="Number in party (7)",
                     "prop_male"="Proportion male (4)",
                     "prop_o65"="Proportion over 65 (1)",
                     "dist_km"="Distance traveled (5)",
                     "ha"="Surface area (8)",
                     "long"="Longitude (10)",
                     "lat"="Latitude (9)",
                     "pop_2010"="County population (10)"))

## Create plot of mean variable importance with 95% CI
## Publication plot
b=ggplot(varimp_plotdat)+
  geom_pointrange(aes(x=mean_imp,xmin=lci95,xmax=uci95,y=name),
                  size=1.5,
                  color="gray30")+
  coord_cartesian(xlim=c(0,0.43),
                  ylim=c(0.4,19.6),
                  expand=F)+
  labs(x="Mean variable importance (95% CI)",
       y="")+
  pubtheme

###############################################################
###############################################################
## Combine variable importance plots
a1=a+
  scale_x_continuous(limits=c(0,0.51),
                     breaks=seq(0,0.5,0.1),
                     name="")+
  annotate("text",label="A",x=0.5,y=1,hjust=1,vjust=0,size=16)+
  coord_cartesian(xlim=c(-0.01,0.51),
                  ylim=c(0,20.6),
                  expand=F)

b1=b+
  scale_x_continuous(breaks=seq(0,0.5,0.1),
                     name="Variable importance")+
  annotate("text",label="B",x=0.5,y=1,hjust=1,vjust=0,size=16)+
  coord_cartesian(xlim=c(-0.01,0.51),
                  ylim=c(0,19.6),
                  expand=F)

## Combine variable importance plots with global on top and species/taxa on bottom
varimp_out=a1/b1

## Export plot
ggsave(plot=varimp_out,"fig3_varimp.png",width=10,height=10,bg="white")

################################################################################
## Create heatmap showing variable importance by species
## Manipulate data
varimp_hm=varimp%>%
  mutate(prop_vi1=case_when(prop_vi==0 ~ NA,
                            TRUE ~ prop_vi),
         name=recode(var,
                     "sppgrp"="Species",
                     "tl"="Length",
                     "target_spp"="Target species",
                     "generalist"="Generalist",
                     "legal"="Legal",
                     "doy"="Day of year",
                     "dayofweek"="Day of week",
                     "day_period"="Time of day",
                     "year"="Year",
                     "bs"="Boat or shore",
                     "hrs"="Hours angling",
                     "instate"="In-state",
                     "tot_ang"="Number in party",
                     "prop_male"="Proportion male",
                     "prop_o65"="Proportion over 65",
                     "dist_km"="Distance traveled",
                     "ha"="Surface area",
                     "long"="Longitude",
                     "lat"="Latitude",
                     "pop_2010"="County population"))%>%
  select(spp,name,prop_vi1)

## Create plot
ggplot(varimp_hm,aes(x=name,y=spp,fill=prop_vi1)) +
  geom_tile(color="white") +
  scale_fill_viridis_c(limits=c(0,1),
                       breaks=seq(0,1,0.2),
                       option="turbo",
                       na.value="gray90",
                       name="Proportional variable importance",
                       guide=guide_colorbar(title.position="top",
                                            title.hjust=0.5,
                                            label.position="bottom",
                                            barwidth=unit(7,"cm"),
                                            barheight=unit(0.7,"cm"),
                                            direction="horizontal"),
                       expand=c(0,0))+
  scale_y_discrete(limits=rev)+
  labs(x="",y="") +
  pubtheme +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        panel.grid=element_blank(),
        legend.position="top",
        legend.box="horizontal",
        legend.title=element_text(size=16),
        legend.text=element_text(size=14))

## Export plot
ggsave(plot=last_plot(),"fig4_propvarimp.png",width=10,height=8,bg="white")

################################################################################
################################################################################
## Summary data for reporting in the manuscript
################################################################################
################################################################################
## Manipulate data to add a unique creel identifier
dat2=dat%>%
  mutate(locyr=paste0(impd,year))

################################################################################
## General summary statistics
## Number of unique creel surveys
length(unique(dat2$locyr))

## Number of unique locations
length(unique(dat2$impd))

## Number of creels per year
tmp=dat2%>%
  group_by(impd,year)%>%
  summarize(n=n())%>%
  ungroup()%>%
  group_by(year)%>%
  summarize(n=n())

min(tmp$n)
max(tmp$n)
mean(tmp$n)
sd(tmp$n)

## Number of creels per impoundment
tmp1=dat2%>%
  group_by(impd,year)%>%
  summarize(n=n())%>%
  ungroup()%>%
  group_by(impd)%>%
  summarize(n=n())

min(tmp1$n)
max(tmp1$n)
mean(tmp1$n)
sd(tmp1$n)

## Number of fish records
nrow(dat2)

################################################################################
## Table 1
## Number of creels each sppgrp was observed
dat2%>%
  group_by(sppgrp,locyr)%>%
  summarize(n=n())%>%
  ungroup()%>%
  group_by(sppgrp)%>%
  summarize(surveys=length(unique(locyr)))%>%
  ungroup()%>%
  mutate(prop=surveys/337)

## Number of each species observed, harvested, and released
dat2%>%
  group_by(sppgrp)%>%
  summarize(n=n(),
            harvest=length(harv[harv=="harvest"]),
            release=n-harvest)

## Bring in data from rpart models to complete table 1
## Training records, testing records, model accuracy
nrow(bcftrain); nrow(bcftest); mean(bcftest$harv==bcftest$pred)
nrow(ccftrain); nrow(ccftest); mean(ccftest$harv==ccftest$pred)
nrow(ccarptrain); nrow(ccarptest); mean(ccarptest$harv==ccarptest$pred)
nrow(crappietrain); nrow(crappietest); mean(crappietest$harv==crappietest$pred)
nrow(fhctrain); nrow(fhctest); mean(fhctest$harv==fhctest$pred)
nrow(fwdtrain); nrow(fwdtest); mean(fwdtest$harv==fwdtest$pred)
nrow(lmbtrain); nrow(lmbtest); mean(lmbtest$harv==lmbtest$pred)
nrow(wipertrain); nrow(wipertest); mean(wipertest$harv==wipertest$pred)
nrow(percidtrain); nrow(percidtest); mean(percidtest$harv==percidtest$pred)
nrow(smbtrain); nrow(smbtest); mean(smbtest$harv==smbtest$pred)
nrow(sunfishtrain); nrow(sunfishtest); mean(sunfishtest$harv==sunfishtest$pred)
nrow(whbtrain); nrow(whbtest); mean(whbtest$harv==whbtest$pred)
nrow(train); nrow(test); mean(test$harv==test$pred)

################################################################################
## Table 2
## Mean and SD for continuous variables
mean(dat1$tl)
sd(dat1$tl)

mean(dat1$doy)
sd(dat1$doy)

mean(dat1$year)
sd(dat1$year)

mean(dat1$hrs,na.rm=T)
sd(dat1$hrs,na.rm=T)

mean(dat1$tot_ang,na.rm=T)
sd(dat1$tot_ang,na.rm=T)

mean(dat1$prop_male,na.rm=T)
sd(dat1$prop_male,na.rm=T)

mean(dat1$prop_o65,na.rm=T)
sd(dat1$prop_o65,na.rm=T)

mean(dat1$dist_km,na.rm=T)
sd(dat1$dist_km,na.rm=T)

mean(dat1$ha,na.rm=T)
sd(dat1$ha,na.rm=T)

mean(dat1$pop_2010,na.rm=T)
sd(dat1$pop_2010,na.rm=T)

mean(dat1$long,na.rm=T)
sd(dat1$long,na.rm=T)

mean(dat1$lat,na.rm=T)
sd(dat1$lat,na.rm=T)