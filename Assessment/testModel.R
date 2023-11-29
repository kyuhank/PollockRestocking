

InputData=InputData[[2]]

chains=10
InputPriorDraws=10^3
#JuvCollapse=52,
#AdultCollapse=62,
nupdate=999
R0_min=2e+6
R0_max=2e+6
Mmedian=0.22                          ## from Kooka (2012)
sigM=0.0
SamNumLimit=300
Steep_alpha=100
Steep_beta=100
gamma_min=1
gamma_max=1
verbose=1
L95nestDiff=0
L05nestDiff=0




SimModel <- cmdstanr::cmdstan_model('../Model/Main/main.stan',
                                    pedantic=F, 
                                    force_recompile=F)

SimObj=SimModel


SamNum=0

sampleStore=list()

InputData$sigR=0

## recruitment deviation sd
sigR=InputData$sigR


sigL05=InputData$sigL05

L05_min=InputData$L05_low
L05_max=InputData$L05_low

sigL95=InputData$sigL95

L95_min=InputData$L95_low
L95_max=InputData$L95_low

JuvCollapse=InputData$JuvCollapseYear
AdultCollapse=InputData$AduCollapseYear

ThresCollapse=InputData$ThresCollapse

nyears=InputData$nyears
A=InputData$A


if(InputData$gamma>0) {
  gamma_min=InputData$gamma
  gamma_max=InputData$gamma
}



RdevCor=InputData$RdevCor

L05devCor=InputData$L05devCor

L95devCor=InputData$L95devCor


Sig=MvnSigmaAR1(sig = sigR, rho = RdevCor, dim = nyears+A-1)

SigL05=MvnSigmaAR1(sig = sigL05, rho = L05devCor, dim = nyears)

SigL95=MvnSigmaAR1(sig = sigL95, rho = L95devCor, dim = nyears)


  ##############################################
  ### Random draws from initial input priors ###
  ##############################################
  
  ## random draws for every iteration from the joint prior (original BH parameterisation)
  InputPrior=draws_array("R0"=replicate(chains, runif(InputPriorDraws, min = R0_min, max = R0_max)),
                         #"R0"=replicate(chains, rlnorm(InputPriorDraws, log(1e+9), 2)),
                         "Natural_M"=replicate(chains, exp(rnorm(InputPriorDraws, log(Mmedian), sigM))),
                         "steepness"=replicate(chains, rbeta(InputPriorDraws, Steep_alpha, Steep_beta)*0.8+0.2),
                         "gamma"=replicate(chains, runif(InputPriorDraws, gamma_min, gamma_max)),
                         "L95"=replicate(chains, runif(InputPriorDraws, L95_min+L95nestDiff, L95_max-L95nestDiff)),
                         "L05"=replicate(chains, runif(InputPriorDraws, L05_min+L05nestDiff, L05_max-L05nestDiff)),
                         .nchains = chains)
  
  
  ### attach Rdev values ####
  #Rdev=aperm( replicate(chains, MASS::mvrnorm(InputPriorDraws, mu = rep(0, dim(Sigma)[1]), Sigma = Sigma), c(2,3,1) ))
  #Rdev=aperm( replicate(chains, replicate(InputPriorDraws, AutoRdev(RdevCor, sigR, nyears+A-1) )), c(2,3,1) )
  
  Rdev=aperm(  replicate(chains, replicate(InputPriorDraws, MASS::mvrnorm(1, mu = rep(0, dim(Sig)[1]), Sigma = Sig))) , c(2,3,1) )
  L05dev=aperm(  replicate(chains, replicate(InputPriorDraws, MASS::mvrnorm(1, mu = rep(0, dim(SigL05)[1]), Sigma = SigL05))) , c(2,3,1) )
  L95dev=aperm(  replicate(chains, replicate(InputPriorDraws, MASS::mvrnorm(1, mu = rep(0, dim(SigL95)[1]), Sigma = SigL95))) , c(2,3,1) )
  
  
  #browser()
  
  dimnames(Rdev)[[3]]<-paste("Rdev[",1:(nyears+A-1),"]", sep='')
  InputPrior=as_draws(abind::abind(InputPrior, Rdev))
  
  dimnames(L05dev)[[3]]<-paste("L05dev[",1:nyears,"]", sep='')
  InputPrior=as_draws(abind::abind(InputPrior, L05dev))
  
  dimnames(L95dev)[[3]]<-paste("L95dev[",1:nyears,"]", sep='')
  InputPrior=as_draws(abind::abind(InputPrior, L95dev))
  
  
  
  #ay=rep(100, 76)
  
  #jy=rep(100, 76)
  
  
  ay=rep(5, 76)
  
  jy=rep(5, 76)
  
  
  InputData$AYt=ay
  
  InputData$JYt=jy
  
  InputData$sigR=0
  
  InputData$RdevCor=0
  
  ############################
  ### Conditioning begins ####
  ############################
  PushForwardCheck<-SimObj$generate_quantities(data = InputData, 
                                                               fitted_params =InputPrior, 
                                                               parallel_chains = dim(InputPrior)[2])
  
  
  
  
  
  
  boxplot(as_draws_matrix(PushForwardCheck$draws("SSBStatus")), ylim=c(0, 2))
  
  
  
  
  
  
