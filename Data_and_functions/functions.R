
## find sigma based on CV in a lognormal distribution
lnSigmaFromCV=function(CV) {
  sqrt(log(CV^2 + 1))
}


## AR1 recruitment deviations
AutoRdev<-function(rho, sigR, times) {
  
  Rdev=c()
  
  Rdev[1]=rnorm(1, 0, sigR)
  
  for (i in 2:times) {
    
    Rdev[i]=rho*Rdev[i-1]+rnorm(1, 0, sigR)
    
  }
  
  Rdev=Rdev-0.5*sigR^2
  
  
  return(Rdev)
  
}


## make Sigma for AR1
MvnSigmaAR1<-function (sig, rho, dim) {
  
  SigVec=rep(sig, dim)
  CorMat=diag(1, dim)
  
  for (i in 1:dim) {
    for (j in 1:dim) {
      CorMat[i,j]=rho^(abs(i-j))
    }
  }
  
  VarCov=diag(SigVec)%*%CorMat%*%diag(SigVec)
  
  return (VarCov) 
  
}


## Selectivity function
Sel=function(bins,
              L05, 
              L95) {
  
  Selec=c()
  
  nbins=length(bins)
  
  for(i in 1:nbins) {
    Selec[i]=19/(19+exp(log(361)*( (bins[i]-L95)/(L05-L95) ) ))
  }
  
  return(Selec)         
}



## Stock recruitment
SRfunction=function(R0,
                    steepness, 
                    SSB, 
                    SSB0, 
                    gamma, 
                    steepPoint=0.2, 
                    Type) {
  
  Rec=c()
  
  for (i in 1:length(SSB)) {
    
    
    if (Type==1) {
      
      alpha=(steepness*R0*(1-steepPoint^(gamma)))/(SSB0^(gamma)*steepPoint^(gamma)*(1-steepness))
      beta=(1/(steepPoint^(gamma))*steepness-1)/(SSB0^(gamma)*(1-steepness))                  
      Rec[i]=(alpha*SSB[i]^(gamma))/(1+beta*SSB[i]^(gamma))
      
    } else if (Type==2) {
      
      if (SSB[i]<SSB0*steepPoint) {
        
        alpha=(steepness*R0*(1-steepPoint^(gamma)))/(SSB0^(gamma)*steepPoint^(gamma)*(1-steepness))
        beta=(1/(steepPoint^(gamma))*steepness-1)/(SSB0^(gamma)*(1-steepness))                  
        Rec[i]=(alpha*SSB[i]^(gamma))/(1+beta*SSB[i]^(gamma))
        
      } else {
        
        alpha=(steepness*R0*(1-steepPoint))/(SSB0*steepPoint*(1-steepness))
        beta=(1/steepPoint*steepness-1)/(SSB0*(1-steepness))                  
        Rec[i]=(alpha*SSB[i])/(1+beta*SSB[i])
        
      }
      
    }
    
    
  }
  
  return(Rec)
  
}


## survival rates at equilibrium
Surv=function(Natural_M,
              A) {
  
  SurAtAge=c()
  
  for (a in 1:A) {
    
    if(a==1) {
      SurAtAge[a]=1
    } else if (a>1 || a<A) {
      SurAtAge[a]=SurAtAge[a-1]*exp(-Natural_M)
    } else if (a==A){
      SurAtAge[a]=SurAtAge[a-1]*(exp(-Natural_M)/(1-exp(-Natural_M)))
    }
    
    
  }
  
  return(SurAtAge) 
}


## egg production per recruit at equilibrium 
eqE=function(A,
             SurAtAge,
             LengthDist,
             Fecundity,
             Maturity,
             FemaleProp) {
  
  eqEa=c()
  
  for (a in 1:A) {
    eqEa[a]=sum(SurAtAge[a]*LengthDist[a,]*Fecundity*Maturity*FemaleProp)
  }
  
  
  return(sum(eqEa))
  
}


############################################################################################################################################################
############################################################################################################################################################

#################################
##### PriorPushForwardCheck #####
#################################

PriorPushForwardCheck=function(InputData,
                               chains=10,
                               InputPriorDraws=10^3,
                               #JuvCollapse=52,
                               #AdultCollapse=62,
                               nupdate=999999,
                               SimObj,
                               R0_min=2e+6,
                               R0_max=2e+10,
                               Mmedian=0.22,                          ## from Kooka (2012)
                               sigM=0.5,
                               SamNumLimit=300,
                               Steep_alpha=1,
                               Steep_beta=1,
                               gamma_min=1,
                               gamma_max=3,
                               verbose=1,
                               L95nestDiff=0,
                               L05nestDiff=0,
                               ...) {
  
  SamNum=0
  
  sampleStore=list()
  
  ## recruitment deviation sd
  sigR=InputData$sigR
  
  
  sigL05=InputData$sigL05
  
  L05_min=InputData$L05_low
  L05_max=InputData$L05_up
  
  sigL95=InputData$sigL95
  
  L95_min=InputData$L95_low
  L95_max=InputData$L95_up
  
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
  
  
  for (i in 1:nupdate) {
    
    ##############################################
    ### Random draws from initial input priors ###
    ##############################################
    
    ## random draws for every iteration from the joint prior (original BH parameterisation)
    InputPrior=draws_array("R0"=replicate(chains, exp(runif(InputPriorDraws, min = log(R0_min), max = log(R0_max)))),
      #"R0"=replicate(chains, rlnorm(InputPriorDraws, log(1e+9), 2)),
      "Natural_M"=replicate(chains, exp(rnorm(InputPriorDraws, log(Mmedian), sigM))),
                           "steepness"=replicate(chains, ifelse(rep(InputData$steepness>0, InputPriorDraws), rep(InputData$steepness, InputPriorDraws), rbeta(InputPriorDraws, Steep_alpha, Steep_beta)*0.8+0.2)),
                           "gamma"=replicate(chains, runif(InputPriorDraws, gamma_min, gamma_max)),
                           "L95"=replicate(chains, runif(InputPriorDraws, L95_min+L95nestDiff, L95_max-L95nestDiff)),
                           "L05"=replicate(chains, runif(InputPriorDraws, L05_min+L05nestDiff, L05_max-L05nestDiff)),
                           .nchains = chains)

    
    ### attach Rdev values ####
<<<<<<< HEAD:Data_and_functions/functions.R
    #Rdev=aperm( replicate(chains, MASS::mvrnorm(InputPriorDraws, mu = rep(0, dim(Sigma)[1]), Sigma = Sigma), c(2,3,1) ))
    #Rdev=aperm( replicate(chains, replicate(InputPriorDraws, AutoRdev(RdevCor, sigR, nyears+A-1) )), c(2,3,1) )
=======
    Rdev=aperm(  replicate(chains, replicate(InputPriorDraws, MASS::mvrnorm(1, mu = rep(-0.5*sigR^2, dim(Sig)[1]), Sigma = Sig))) , c(2,3,1) )
    L05dev=aperm(  replicate(chains, replicate(InputPriorDraws, MASS::mvrnorm(1, mu = rep(0, dim(SigL05)[1]), Sigma = SigL05))) , c(2,3,1) )
    L95dev=aperm(  replicate(chains, replicate(InputPriorDraws, MASS::mvrnorm(1, mu = rep(0, dim(SigL95)[1]), Sigma = SigL95))) , c(2,3,1) )
>>>>>>> 75624bdab81366e7ef76bd5316ef91d3e94956f2:data/functions.R
    
    Rdev=aperm(  replicate(chains, replicate(InputPriorDraws, MASS::mvrnorm(1, mu = rep(-0.5*sigR^2, dim(Sig)[1]), Sigma = Sig))) , c(2,3,1) )
    L05dev=aperm(  replicate(chains, replicate(InputPriorDraws, MASS::mvrnorm(1, mu = rep(-0.5*sigL05^2, dim(SigL05)[1]), Sigma = SigL05))) , c(2,3,1) )
    L95dev=aperm(  replicate(chains, replicate(InputPriorDraws, MASS::mvrnorm(1, mu = rep(-0.5*sigL95^2, dim(SigL95)[1]), Sigma = SigL95))) , c(2,3,1) )
    
    
    #browser()
    
    dimnames(Rdev)[[3]]<-paste("Rdev[",1:(nyears+A-1),"]", sep='')
    InputPrior=as_draws(abind::abind(InputPrior, Rdev))
    
    dimnames(L05dev)[[3]]<-paste("L05dev[",1:nyears,"]", sep='')
    InputPrior=as_draws(abind::abind(InputPrior, L05dev))
    
    dimnames(L95dev)[[3]]<-paste("L95dev[",1:nyears,"]", sep='')
    InputPrior=as_draws(abind::abind(InputPrior, L95dev))
    
      
    ############################
    ### Conditioning begins ####
    ############################
    capture.output({PushForwardCheck<-SimObj$generate_quantities(data = InputData, 
                                                                 fitted_params =InputPrior, 
                                                                 parallel_chains = dim(InputPrior)[2], ...)}, file=nullfile())
    
    
    #PushForwardCheck<-SimObj$generate_quantities(data = InputData, 
    #                                                             fitted_params =InputPrior, 
    #                                                             parallel_chains = dim(InputPrior)[2], ...)
    
    
    
    
    ## draw samples of interest for conditioning
    JBStatus=PushForwardCheck$draws("JBStatus")
    ABStatus=PushForwardCheck$draws("ABStatus")
    BStatus=PushForwardCheck$draws("TBStatus")
    SSBStatus=PushForwardCheck$draws("SSBStatus")
    
   # browser()
    
  #  print(as_draws_matrix(BStatus )[,10])
    
    ## fishery collapse defined based on the assumption: B/B0<0.1 (Worm and Hilborn, 2009)
    SelectedSamples=which(c(JBStatus[,,JuvCollapse])<ThresCollapse & c(ABStatus[,,AdultCollapse])<ThresCollapse 
                          & apply(as_draws_matrix(JBStatus[,,(1:nyears)]), 1, function(x) all(x>=0.0)) 
                          & apply(as_draws_matrix(ABStatus[,,(1:nyears)]), 1, function(x) all(x>=0.0))
                          & apply(as_draws_matrix(SSBStatus[,,(1:nyears)]), 1, function(x) all(x<=1.5)) )
    
    ## filtering the samples
    PriorSample=as_draws_matrix(InputPrior)[SelectedSamples,]
    
    SamNum=SamNum+length(SelectedSamples)
    
    sampleStore[[i]]=PriorSample
    
    ## report the process
    if(verbose==1) { 
      print(paste("iter-",i,": ", SamNum, "/", SamNumLimit, sep=""))
      #print(paste("run-",i,": ", length(SelectedSamples), "/", InputPriorDraws*chains, " (total: ", SamNum,")", sep=""))
    }
    
    ## stopping criteria  
     if(SamNum>=SamNumLimit | i==nupdate ) {
       sampleStore=as_draws_array(do.call(rbind, sampleStore))
       sampleStore=sampleStore[1:SamNumLimit,,]
       print(paste(SamNum," samples obtained!", sep=""))
      
       break
      
     }
    
    
  }
  
  return("UpdatedSamples"=as_draws_df(sampleStore))
  
  }


###############################
##### Projection function #####
###############################

DeriveQuantities=function(InputData,
                          FittedParams,
                          SimObj,
                          nRelease,
                          nReleasePeriod=10,
                          ReleaseStart=70) {
  
  
  nProj=InputData$nProj
  nCatchYears=InputData$nCatchYears
  
  Release=rep(0, nCatchYears+nProj)
  
  if(nRelease>0) {
    Release[ReleaseStart:(ReleaseStart+nReleasePeriod-1)]=rep(nRelease,nReleasePeriod)
  }
  
  Inputs=InputData
  
  Inputs$Release=Release
  
  capture.output({ModelPredicted<-SimObj$generate_quantities(data=Inputs,
                                                             fitted_params = FittedParams)}, file=nullfile())
  
  
  
  #ABt=as_draws_matrix(ModelPredicted$draws("ABt"))
  #JBt=as_draws_matrix(ModelPredicted$draws("JBt"))
  TJBt=as_draws_matrix(ModelPredicted$draws("TJBt"))
  
  #Bt=as_draws_matrix(ModelPredicted$draws("Bt"))
  #availBt=as_draws_matrix(ModelPredicted$draws("availBt"))
  
  SSBrel=as_draws_matrix(ModelPredicted$draws("SSBrel"))
  JBrel=as_draws_matrix(ModelPredicted$draws("JBrel"))
  ABrel=as_draws_matrix(ModelPredicted$draws("ABrel"))
  
  
  
  #JBStatus=as_draws_matrix(ModelPredicted$draws("JBStatus"))
  #ABStatus=as_draws_matrix(ModelPredicted$draws("ABStatus"))
#  availBStatus=as_draws_matrix(ModelPredicted$draws("availBStatus"))
  #TBStatus=as_draws_matrix(ModelPredicted$draws("TBStatus"))
  
#  JHt=as_draws_matrix(ModelPredicted$draws("JHt"))
#  AHt=as_draws_matrix(ModelPredicted$draws("AHt"))
  
#  JCt=as_draws_matrix(ModelPredicted$draws("JCt"))
#  ACt=as_draws_matrix(ModelPredicted$draws("ACt"))
  
  SSBt=as_draws_matrix(ModelPredicted$draws("SSBt"))
  SSBStatus=as_draws_matrix(ModelPredicted$draws("SSBStatus"))
  
  sf <- ModelPredicted$summary(variables = NULL,'mean', ~quantile(.x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm=T),'sd')
  
  return(list("sf"=sf,
              #"ABt"=ABt,
              "TJBt"=TJBt,
              #"Bt"=Bt,
              "SSBrel"=SSBrel,
              "JBrel"=JBrel,
              "ABrel"=ABrel,
              #"availBt"=availBt,
              #"JBStatus"=JBStatus,
              #"ABStatus"=ABStatus,
              "SSBStatus"=SSBStatus,
#              "availBStatus"=availBStatus,
              #"TBStatus"=TBStatus,
#              "JHt"=JHt,
#              "AHt"=AHt,
#              "JCt"=JCt,
#              "ACt"=ACt,
              "SSBt"=SSBt) )
  
}
