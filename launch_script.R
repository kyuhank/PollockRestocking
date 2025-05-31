#--------------------------------------------------------------
# Length-Based Age-Structured Model by Kyuhan Kim
# Copyright © 2023 Kyuhan Kim. All rights reserved.
# Contact: kh2064@gmail.com for questions
# MIT License: https://opensource.org/licenses/MIT
# --------------------------------------------------------------

# This script is to run the pollock assessment
# run this entire script to run the job
# it will take days to finish the job, so it was originally designed to run on gateaux.io (api server)

##############
## Preamble ##
##############

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

require(bakeR)
require(dplyr)

## only local run mode works with this version (original version was designed for gateaux.io)
LocalRun=T 

#####################
## sampling option ##
#####################

nDraws=400
nChains=30
minSampleNum=200

#######################################
## Input values for making scenarios ##
#######################################

## release options ##
## base case: average annual release number (230,000 fish) by South Korea gov
baseRelease=230000
multiple=c(0, 1, 10, 100, 1000, 10000)
nRelease=multiple*baseRelease

## input pars for making scenarios (trigger parallel jobs) ##
gamma=c(1,2)                ## can also put prior (e.g., gamma=-99. gamma_min; gamma_max=3)
steepness=c(-99)              ## can also put prior (e.g., steep=-99. Steep_alpha; Steep_beta)

RdevCor=c(0, 0.4)
L05devCor=c(0.0)              # no autocorrelation
L95devCor=c(0.0)              # no autocorrelation


JuvCollapseYear=c(52, 55)         ## 52: 1997;  55: 2000;  58: 2003;
AduCollapseYear=c(55, 58)

PmixL=c(0.5, 5)

## nested runs (run within each job)


JuvMaxVul=c(0.5, 1)
AduMaxVul=c(0.5, 1)
sigR=1                      ## from US assessment reports
sigR_cond=sigR * sqrt(1-RdevCor^2)  

sigL05=c(0.0)               # no autocorrelation
sigL95=c(0.0)               # no autocorrelation
NatExponent=c(0.0, -1)      ## Lorenzen et al., 2020

## prior setting ##
R0_min=c(2e+8)             ## from Kang et al., 2013
R0_max=c(2e+11)            ## 1000 times larger than R0_min
Steep_alpha=2          
Steep_beta=2

## fix these values ##
L05_low=10                 ## from Kang et al., 2013
L05_up=15                  ## from Kang et al., 2013
L95_low=20                 ## from Kang et al., 2013
L95_up=25                  ## from Kang et al., 2013
Mmedian=0.22               ## from Kooka et al., 2012
sigM=0.3                 
gamma_min=1
gamma_max=3

########################
#### pre-processing ####
########################

#### all possible combinations ####
parsComb=merge(expand.grid("gamma"=gamma,
                           "steepness"=steepness,
                           "RdevCor"=RdevCor,
                           "L05devCor"=L05devCor,
                           "L95devCor"=L95devCor,
                           "JuvCollapseYear"=JuvCollapseYear,
                           "AduCollapseYear"=AduCollapseYear,
                           "PmixL"=PmixL),
               cbind(R0_min, R0_max, Steep_alpha, Steep_beta, L05_low, L05_up, L95_low, L95_up, gamma_min, gamma_max, Mmedian, sigM))


parsSensitivity <- lapply(1:nrow(parsComb), function(ll){
  list(pars = c(as.list(parsComb[ll,]), as.list(c("nDraws"=nDraws,
                                                  "nChains"=nChains,
                                                  "minSampleNum"=minSampleNum,
                                                  "NatExponent"=NatExponent,
                                                  "sigR"=sigR,
                                                  "sigL05"=sigL05,
                                                  "sigL95"=sigL95,
                                                  "JuvMaxVul"=JuvMaxVul,
                                                  "AduMaxVul"=AduMaxVul
                                                  #"AdultFraction"=AdultFraction
  )),
  nRelease=as.list(nRelease), ReleaseScenarioNumber=length(nRelease)) )
})

names(parsSensitivity) <- paste('RUN',1:length(parsSensitivity), sep='_')

##########################################################################################
################################### run the job (SCRA) ###############################
##########################################################################################

if(LocalRun!=1) {
  
  if (!exists("JWT")) {
    stop("API token is not provided")
  } else {
    
    gateaux_job_runner(parsSensitivity,
                       server = "gateaux.io/api", 
                       JWT=JWT,
                       report_name = "PollockRelease-analysis",
                       log_jobs = F)
  }
} else {
  source("assessment.R")
}
