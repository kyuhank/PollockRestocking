##############
## Preamble ##
##############

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

require(bakeR)
require(dplyr)

LocalRun=F


#following jobs were rerun because of unnknown error
#reRunIndex=which(c(4628:4659) %in% c(4628, 4629, 4633, 4653))

#####################
## sampling option ##
#####################

nDraws=300
nChains=14
minSampleNum=200

#######################################
## Input values for making scenarios ##
#######################################

## release options ##
## base case: 10 times approx release number (200,000 fish) by South Korea gov
baseRelease=2e+6
multiple=c(0, 1, 100, 200, 300, 400, 500, 600, 700)
nRelease=multiple*baseRelease

## input pars for making scenarios (trigger parallel jobs) ##
gamma=c(1, 2)                ## can also put prior (e.g., gamma=-99; gamma_min=1; gamma_max=3)
RdevCor=c(0, 0.8)
L05devCor=c(0.0, 0.8)
L95devCor=c(0.0)

JuvCollapseYear=c(52, 55)         ## 52: 1997;  55: 2000;  58: 2003;
AduCollapseYear=c(55, 58)         


## nested runs (run within each job)
PmixL=c(0.5, 5)
AdultFraction=c(0.2, 0.8)
JuvMaxVul=c(0.5, 1)
AduMaxVul=c(0.5, 1)
sigR=c(1)                  ## from US assessment reports
sigL05=c(0.0)
sigL95=c(0.0)
NatExponent=c(0.0,-1)


## prior setting ##
R0_min=c(2e+6)
R0_max=c(2e+10)
Steep_alpha=2
Steep_beta=2

## fix these values ##
L05_low=13
L05_up=17
L95_low=23
L95_up=27
Mmedian=0.22
sigM=0.5
gamma_min=1
gamma_max=3

########################
#### pre-processing ####
########################

#### all possible combinations ####
parsComb=merge(expand.grid("gamma"=gamma, 
                           "RdevCor"=RdevCor,
                           "L05devCor"=L05devCor,
                           "L95devCor"=L95devCor,
                           "JuvCollapseYear"=JuvCollapseYear,
                           "AduCollapseYear"=AduCollapseYear),
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
                                             "AduMaxVul"=AduMaxVul,
                                             "AdultFraction"=AdultFraction,
                                             "PmixL"=PmixL
                                             )),
                nRelease=as.list(nRelease), ReleaseScenarioNumber=length(nRelease)) )
})

names(parsSensitivity) <- paste('RUN',1:length(parsSensitivity), sep='_')

parsSensitivity=parsSensitivity[reRunIndex]

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

