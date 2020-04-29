rm(list=ls())

install.packages("writexl")

library(metafor)
library(readxl)
library(tidyverse)
library(writexl)

ds <- read_excel('final data.xlsx')
as.data.frame(table(ds$`Mitochondrial function and iron metabolism measurements`))->dsS
variables <- dsS$Var1
var <- droplevels(variables[-31])
#do something about DRP1 protein expression [var 31]

###### Meta-analysis###############

fits<-list()
for(v in c(1:length(var))){
  rma(m1i=Mean.frail,
      m2i=Mean.ctrl,
      sd1i=SD.frail,
      sd2i=SD.ctrl,
      n1i=N.frail,
      n2i=N.ctrl,
      slab=Study,
      method='REML',
      measure='SMD',
      data=ds[ds$`Mitochondrial function and iron metabolism measurements`==var[v],])->fit
  fits[[v]]<-fit}
names(fits)<-var

#DRP-1 protein expression (variable no. 31) didn't fit with the REML model so we used DL.
rma(m1i=Mean.frail,
    m2i=Mean.ctrl,
    sd1i=SD.frail,
    sd2i=SD.ctrl,
    n1i=N.frail,
    n2i=N.ctrl,
    slab=Study,
    method='DL',
    measure='SMD',
    data=ds[ds$`Mitochondrial function and iron metabolism measurements`==variables[31],])->DRP
###############################

lapply(fits,function(x){
  c('estimate'=x$b,
    'pval'=x$pval,
    'Q.pval'=x$QEp,
    'I^2'=x$I2,
    'upper CI'=x$ci.ub,
    'lower CI'=x$ci.lb)})->fits.pvals
do.call(rbind,fits.pvals)->ds.models
as.data.frame(ds.models)->ds.models

pbl_bias <- list()
for (f in fits){
  pb <- regtest.rma(f)
  pbl_bias <- c(append(pbl_bias, pb$pval))}
names(pbl_bias)<-var
pbl_bias_frame<-as.data.frame(pbl_bias)
pbl_bias_final<-t(pbl_bias_frame)

names <- row.names(ds.models)
Final_data<-cbind(names,pbl_bias_final, ds.models)
write_xlsx(Final_data, "estimates.xlsx")
#dont forget to add DRP1 manually to the data later
################################

##### Moderator analysis ######

#Spiecies as a moderator
species_fits<-list()
rma(m1i=Mean.frail,
    m2i=Mean.ctrl,
    sd1i=SD.frail,
    sd2i=SD.ctrl,
    n1i=N.frail,
    n2i=N.ctrl,
    slab=Study,
    method='REML',
    measure='SMD',
    data=ds[ds$`Mitochondrial function and iron metabolism measurements`==var[80],],
    mods = ~ Species-1)->fit
species_fits<-c(append(species_fits,fit$QMp))
new_var<- droplevels(var[-c(1,6,9,18,22,44,49,51,75,77)])
moderator_test_S <- as.data.frame(species_fits)
moderator_test_S <- t(moderator_test_S)
names<-as.data.frame(new_var)
moderator_test_S<-cbind(names,moderator_test_S)
write_xlsx(moderator_test_S, "moderator species.xlsx")

#Type of frailty assessment as a moderator
frailty_fits<-list()
rma(m1i=Mean.frail,
    m2i=Mean.ctrl,
    sd1i=SD.frail,
    sd2i=SD.ctrl,
    n1i=N.frail,
    n2i=N.ctrl,
    slab=Study,
    method='REML',
    measure='SMD',
    data=ds[ds$`Mitochondrial function and iron metabolism measurements`==var[80],],
    mods = ~ Frailty.assesment.model-1)->fit
frailty_fits<-c(append(frailty_fits,fit$QMp))

other_new_var<- droplevels(var[-c(9,10,18,20,38,41,45,46,48,49,53,67,80)])
moderator_test_F <- as.data.frame(frailty_fits)
moderator_test_F <- t(moderator_test_F)
names<-as.data.frame(other_new_var)
moderator_test_F<-cbind(names,moderator_test_F)
write_xlsx(moderator_test_F, "moderator frailty assessment.xlsx")
###################################

######## Subgroup analysis ########
ds <- read_excel('final data.xlsx')
variable<-ds %>% filter(`Mitochondrial function and iron metabolism measurements`=='DRP1 gene expression')
species_list<-variable %>% select(Species) %>% distinct() %>% pull(Species)

fits<-list()
for(s in species_list){
  fit<-rma(m1i=Mean.frail, m2i=Mean.ctrl, sd1i=SD.frail, sd2i=SD.ctrl, n1i=N.frail, n2i=N.ctrl,
           slab=Study, method='REML', measure='SMD', data=filter(variable,Species==s))
  fits[[s]]<-fit}

fits<-list(fit_mice,fit_rat,fit_human)
lapply(fits,function(x){
  c('estimate'=x$b,
    'pval'=x$pval,
    'Q.pval'=x$QEp,
    'I^2'=x$I2,
    'upper CI'=x$ci.ub,
    'lower CI'=x$ci.lb)})->fits.pvals
do.call(rbind,fits.pvals)->subgroup.models
as.data.frame(subgroup.models)->subgroup.models
names<-as.data.frame(species_list)
final_data<-cbind(names,subgroup.models)


rm(list=ls())