#--------------------------------------------------------------
# Length-Based Age-Structured Model by Kyuhan Kim
# Copyright © 2023 Kyuhan Kim. All rights reserved.
# Contact: kh2064@gmail.com for questions
# MIT License: https://opensource.org/licenses/MIT
# --------------------------------------------------------------

##############
## Preamble ##
##############

library(posterior)
library(reshape2)
library(ggplot2)
library(cowplot)
library(dplyr)
library(GGally)
library(tidyverse)
library(parallel)


#### input values to make scenarios ####
if (!interactive()) {
  
  LocalRun <- as.integer(Sys.getenv("LocalRun", 0))
  
  ## sampling option ##
  nDraws <- as.integer(Sys.getenv("nDraws"))
  nChains <- as.integer(Sys.getenv("nChains"))
  minSampleNum <- as.numeric(Sys.getenv("minSampleNum"))
  
  ## input pars (parallel run: multiple jobs) ##
  gamma <- as.numeric(Sys.getenv("gamma"))
  steepness <- as.numeric(Sys.getenv("steepness"))
  
  RdevCor <- as.numeric(Sys.getenv("RdevCor"))
  L05devCor <- as.numeric(Sys.getenv("L05devCor"))
  L95devCor <- as.numeric(Sys.getenv("L95devCor"))
  L05devCor <- as.numeric(Sys.getenv("L05devCor"))
  L95devCor <- as.numeric(Sys.getenv("L95devCor"))
  R0_min <- as.numeric(Sys.getenv("R0_min"))
  R0_max <- as.numeric(Sys.getenv("R0_max"))
  Steep_alpha <- as.numeric(Sys.getenv("Steep_alpha"))
  Steep_beta <- as.numeric(Sys.getenv("Steep_beta"))
  L05_low <- as.numeric(Sys.getenv("L05_low"))
  L05_up <- as.numeric(Sys.getenv("L05_up"))
  L95_low <- as.numeric(Sys.getenv("L95_low"))
  L95_up <- as.numeric(Sys.getenv("L95_up"))
  gamma_min <- as.numeric(Sys.getenv("gamma_min"))
  gamma_max <- as.numeric(Sys.getenv("gamma_max"))
  
  JuvCollapseYear <- as.numeric(Sys.getenv("JuvCollapseYear"))
  AduCollapseYear <- as.numeric(Sys.getenv("AduCollapseYear"))
  
  ## for nested runs (run within each job) ##
  sigR <- as.numeric(Sys.getenv()[grep(names(Sys.getenv()), pattern = "sigR" )])
  sigL05 <- as.numeric(Sys.getenv()[grep(names(Sys.getenv()), pattern = "sigL05" )])
  sigL95 <- as.numeric(Sys.getenv()[grep(names(Sys.getenv()), pattern = "sigL95" )])
  JuvMaxVul <- as.numeric(Sys.getenv()[grep(names(Sys.getenv()), pattern = "JuvMaxVul" )])
  AduMaxVul <- as.numeric(Sys.getenv()[grep(names(Sys.getenv()), pattern = "AduMaxVul" )])
  
  #AdultFraction <- as.numeric(Sys.getenv()[grep(names(Sys.getenv()), pattern = "AdultFraction" )])
  
  #CatchSepYr <- as.numerir(Sys.getenv()[grep(names(Sys.getenv()), pattern = "CatchSepYr" )])
  
  
  PmixL <- as.numeric(Sys.getenv()[grep(names(Sys.getenv()), pattern = "PmixL" )])
  NatExponent <- as.numeric(Sys.getenv()[grep(names(Sys.getenv()), pattern = "NatExponent" )])
  Mmedian <- as.numeric(Sys.getenv("Mmedian"))
  sigM <- as.numeric(Sys.getenv("sigM"))
  
}

#### input values to make scenarios ####
if (!interactive()) {
  nRelease=c()
  
  ReleaseScenarioNumber<- as.integer(Sys.getenv("ReleaseScenarioNumber"))
  
  for (i in 1:ReleaseScenarioNumber) {
    nRelease[i] <- as.numeric(Sys.getenv(paste("nRelease",i,sep='')))
  }
  
} 


cat('\n\n Setup finished \n\n')  

######################
### load functions ###
######################
source("../Data_and_functions/data.R")
source("../Data_and_functions/functions.R")

#######################
## compile the model ##
#######################

## simulation model
SimModel <- cmdstanr::cmdstan_model('../Model/Main/main.stan',
                                    pedantic=F, 
                                    force_recompile=F)

cat('\n\n Stan model compilation successful \n\n')

##########################################################################################################################################
##########################################################################################################################################

####################
## Make Input Obj ##
####################

if(LocalRun==1) {
  ExpandRuns=merge(expand.grid("JuvMaxVul"=JuvMaxVul,
                               "AduMaxVul"=AduMaxVul,
                               "sigR"=sigR,
                               "sigL05"=sigL05,
                               "sigL95"=sigL95,
                               "gamma"=gamma,
                               "steepness"=steepness,
                               "RdevCor"=RdevCor,
                               "L05devCor"=L05devCor,
                               "L95devCor"=L95devCor,
                               "JuvCollapseYear"=JuvCollapseYear,
                               "AduCollapseYear"=AduCollapseYear,
                               #"AdultFraction"=AdultFraction,
                               #"CatchSepYr"=CatchSepYr,
                               "PmixL"=PmixL,
                               "NatExponent"=NatExponent),
                   ## not expanding these below
                   cbind("L05_low"=L05_low,
                         "L05_up"=L05_up,
                         "L95_low"=L95_low,
                         "L95_up"=L95_up))
} else {
  ExpandRuns=merge(expand.grid("JuvMaxVul"=JuvMaxVul,
                               "AduMaxVul"=AduMaxVul,
                               "sigR"=sigR,
                               "sigL05"=sigL05,
                               "sigL95"=sigL95,
                               #"CatchSepYr"=CatchSepYr,
                               #"AdultFraction"=AdultFraction,
                               "PmixL"=PmixL,
                               "NatExponent"=NatExponent),
                   ## not expanding these below
                   cbind("gamma"=gamma,
                         "steepness"=steepness,
                         "RdevCor"=RdevCor,
                         "L05devCor"=L05devCor,
                         "L95devCor"=L95devCor,
                         "L05_low"=L05_low,
                         "L05_up"=L05_up,
                         "L95_low"=L95_low,
                         "L95_up"=L95_up,
                         "Mmedian"=Mmedian,
                         "sigM"=sigM,
                         "JuvCollapseYear"=JuvCollapseYear,
                         "AduCollapseYear"=AduCollapseYear))
}


InputData=mapply(MakeInputObj, 
                 sigR=ExpandRuns[,"sigR"],
                 sigL05=ExpandRuns[,"sigL05"],
                 sigL95=ExpandRuns[,"sigL95"],
                 JuvMaxVul = ExpandRuns[,"JuvMaxVul"],
                 AduMaxVul = ExpandRuns[,"AduMaxVul"],
                 L05_low = ExpandRuns[,"L05_low"],
                 L05_up = ExpandRuns[,"L05_up"],
                 L95_low = ExpandRuns[,"L95_low"],
                 L95_up = ExpandRuns[,"L95_up"],
                 gamma=ExpandRuns[,"gamma"],
                 steepness=ExpandRuns[,"steepness"],
                 RdevCor=ExpandRuns[,"RdevCor"],
                 L05devCor=ExpandRuns[,"L05devCor"],
                 L95devCor=ExpandRuns[,"L95devCor"],
                 #AdultFraction=ExpandRuns[,"AdultFraction"],
                 #CatchSepYr=ExpandRuns[,"CatchSepYr"],
                 PmixL=ExpandRuns[,"PmixL"],
                 NatExponent=ExpandRuns[,"NatExponent"],
                 JuvCollapseYear=ExpandRuns[,"JuvCollapseYear"],
                 AduCollapseYear=ExpandRuns[,"AduCollapseYear"],
                 SIMPLIFY = F)

#############################
#### rejection sampling  ####
#############################

cat('\n\n run SRA using rejection sampling \n\n')  

#### Find updated priors ####

UpdatedInputSample=list()

print(paste("nRuns: ",length(InputData),sep=""))

for (i in 1:length(InputData)) {
  
  print(paste("Run-",i,sep=""))
  
  UpdatedInputSample[[i]]=PriorPushForwardCheck(InputData=InputData[[i]],
                                                SimObj=SimModel,
                                                chains=nChains,
                                                InputPriorDraws=nDraws,
                                                SamNumLimit = minSampleNum,
                                                R0_min=R0_min,
                                                R0_max=R0_max,
                                                Steep_alpha=Steep_alpha,
                                                Steep_beta=Steep_beta,
                                                gamma_min = gamma_min,
                                                gamma_max = gamma_max,
                                                Mmedian=Mmedian,
                                                sigM=sigM,
                                                verbose=T)
}


#####################################################
#### projection under assumed release scenarios  ####
#####################################################

cat('\n\n obtain derived quantities \n\n')  

sf=list()
#ABt=list()
TJBt=list()
#Bt=list()
#JBStatus=list()
#ABStatus=list()
#TBStatus=list()
SSBStatus=list()
#availBStatus=list()
#JHt=list()
#AHt=list()
#JCt=list()
#ACt=list()
SSBt=list()
#availBt=list()
SSBrel=list()
JBrel=list()
ABrel=list()

print(paste("nRuns: ",length(InputData),sep=""))

for (i in 1:length(InputData)) {
  
  print(paste("Run-",i,sep=""))
  ## parallel run for each release scenarios
  temp=parallel::mcmapply(DeriveQuantities,
                          mc.cores = length(nRelease),
                          nRelease=nRelease,
                          MoreArgs = list(SimObj = SimModel,
                                          FittedParams = as_draws_array(UpdatedInputSample[[i]]),
                                          InputData= InputData[[i]]),
                          SIMPLIFY = T)
  
  BioQuantOfInterest=c("SSBStatus", "SSBt", "TJBT", "SSBrel", "JBrel", "ABrel")
  
  sf[[i]]=temp["sf",] %>% lapply(function(x) x %>% filter(grepl( paste(BioQuantOfInterest, collapse="|"), variable) ) )
  #TABt[[i]]=temp["TABt",]
  TJBt[[i]]=temp["TJBt",]
  #Bt[[i]]=temp["Bt",]
  SSBrel[[i]]=temp["SSBrel",]
  JBrel[[i]]=temp["JBrel",]
  ABrel[[i]]=temp["ABrel",]
  
  #  availBt[[i]]=temp["availBt",]
  #TBStatus[[i]]=temp["TBStatus",]
  #JBStatus[[i]]=temp["JBStatus",]
  #ABStatus[[i]]=temp["ABStatus",]
  SSBStatus[[i]]=temp["SSBStatus",]
  #  availBStatus[[i]]=temp["availBStatus",]
  
  #  JHt[[i]]=temp["JHt",]
  #  AHt[[i]]=temp["AHt",]
  #  JCt[[i]]=temp["JCt",]
  #  ACt[[i]]=temp["ACt",]
  SSBt[[i]]=temp["SSBt",]
  
}

rm(temp)

##########################################################################################################
##########################################################################################################

###################
##### grooming ####
###################

inputframe <- data.frame("sigR"=ExpandRuns[,"sigR"],
                         "sigL05"=ExpandRuns[,"sigL05"], 
                         "sigL95"=ExpandRuns[,"sigL95"], 
                         "JuvMaxVul"=ExpandRuns[,"JuvMaxVul"], 
                         "AduMaxVul"=ExpandRuns[,"AduMaxVul"],
                         "gamma"=ExpandRuns[,"gamma"],
                         "steepness"=ExpandRuns[,"steepness"],
                         "RdevCor"=ExpandRuns[,"RdevCor"], 
                         "L05devCor"=ExpandRuns[,"L05devCor"],
                         "L95devCor"=ExpandRuns[,"L95devCor"], 
                         "JuvCollapseYear"=ExpandRuns[,"JuvCollapseYear"],
                         "AduCollapseYear"=ExpandRuns[,"AduCollapseYear"],
                         "PmixL"=ExpandRuns[,"PmixL"], 
                         #"CatchSepYr"=ExpandRuns[,"CatchSepYr"],
                         #"AdultFraction"=ExpandRuns[,"AdultFraction"], 
                         "R0_min"=R0_min, 
                         "R0_max"=R0_max, 
                         "Mmedian"=Mmedian,
                         "NatExponent"=ExpandRuns[,"NatExponent"],
                         "sigM"=sigM,
                         "Steep_alpha"=Steep_alpha, 
                         "Steep_beta"=Steep_beta, 
                         "L05_low"=ExpandRuns[,"L05_low"], 
                         "L05_up"=ExpandRuns[,"L05_up"], 
                         "L95_low"=ExpandRuns[,"L95_low"], 
                         "L95_up"=ExpandRuns[,"L95_up"], 
                         "gamma_min"=gamma_min, 
                         "gamma_max"=gamma_max)


for (i in 1:length(InputData)) {
  sf[[i]]=sf[[i]] %>% bind_rows(.id = "ReleaseScenario" ) %>% bind_cols(inputframe[i,])
}


InputData=InputData[[1]]

if (!interactive()) {
  save.image("ModelResults.RData")
} else {
  save.image("assmt_output/ModelResults.RData")
}
