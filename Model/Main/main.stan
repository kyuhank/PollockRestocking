functions{
  #include "../Submodel/Selex.stan"
  #include "../Submodel/StockRecruitment.stan"
  #include "../Submodel/eqQuant.stan"
  #include "../Submodel/dynamics.stan"
  #include "../Submodel/HelperFunctions.stan"
  }

data {
  int nbins;
  int nCatchYears;
  int A;
  int nProj;
  int nyears;
  int nplatoons;
  int FindBrps;
  
  vector[nbins] bins;
  vector[nCatchYears] JYt;
  vector[nCatchYears] AYt;
  vector[nbins] LengthWeight;
  vector[nbins] Maturity;
  vector[nyears] Release;
  vector[nplatoons] Rprop;
  
  real FemaleProp;
  
  real LrefForM;
  
  array[nplatoons] matrix[A, nbins] LengthDist;
  
  real steepPoint;
  //real gamma;
  
  real NatExponent;
  
  //real L05;
  //real L95;
  real sigR;
  real JuvMaxVul;
  real AduMaxVul;
  
  real L05_low;
  real L05_up;
  
  real L95_low;
  real L95_up;
  
  }

parameters {
  
  real R0;
  real Natural_M;
  real steepness;
  real gamma;
  real L95;
  real L05;
  vector[nyears+A-1] Rdev;
  vector[nyears] L05dev;
  vector[nyears] L95dev;
  
  }
  
transformed parameters {
  
  real L05_logit=bounded_logit(L05, L05_low, L05_up);
  
  vector[nyears] L05t_logit=L05_logit+L05dev;
  
  
  real L95_logit=bounded_logit(L95, L95_low, L95_up);
  
  vector[nyears] L95t_logit=L95_logit+L95dev;
  
  vector[nbins] Nat_i;
  
  for (i in 1:nbins) {
    Nat_i[i]=Natural_M*(bins[i]/LrefForM)^NatExponent;
  }  
    
}

generated quantities {
  
  vector[nyears] L05t;
  vector[nyears] L95t;
  matrix[nbins, nyears] JVul;
  vector[nbins] AVul=Maturity*AduMaxVul;
  
  vector[nbins] JVulMean=Sel(bins, L05, L95, nbins).*((1-Maturity)*JuvMaxVul);
  
  for (t in 1:nyears){
    L05t[t]=bounded_inv_logit(L05t_logit[t],  L05_low, L05_up);
    L95t[t]=bounded_inv_logit(L95t_logit[t],  L95_low, L95_up);
    JVul[,t]=Sel(bins, L05t[t], L95t[t], nbins).*((1-Maturity)*JuvMaxVul);
  }
  
  
  matrix[A+16, nyears] OutPut=DynamicsN(nCatchYears,
                                        A,
                                        nplatoons,
                                        bins,
                                        LengthWeight,
                                        Nat_i,
                                        Maturity,
                                        JVulMean,
                                        JVul,
                                        AVul,
                                        R0,
                                        sigR,
                                        steepness,
                                        gamma,
                                        FemaleProp,
                                        LengthDist,
                                        JYt,
                                        AYt,
                                        steepPoint,
                                        nbins,
                                        Release,
                                        Rdev,
                                        Rprop,
                                        nProj
                                        );
  
  matrix[A, nyears] Nat=OutPut[1:A,];
  vector[nyears] JBt=OutPut[A+1,]';
  vector[nyears] ABt=OutPut[A+2,]';
  vector[nyears] SSBt=OutPut[A+3,]';
  vector[nyears] Bt=OutPut[A+4,]';
  vector[nyears] TBStatus=OutPut[A+5,]';
  vector[nyears] JBStatus=OutPut[A+6,]';
  vector[nyears] ABStatus=OutPut[A+7,]';
  vector[nyears] JHt=OutPut[A+8,]';
  vector[nyears] AHt=OutPut[A+9,]';
  vector[nyears] availBt=OutPut[A+10,]';
  
  real SSB0_sim=OutPut[A+11,1];
  
  vector[nyears] ReleaseProp=OutPut[A+12,]';
  
  //real relRelease=max(Release)/R0;
  
  
  vector[nyears] JCt=OutPut[A+13,]';
  vector[nyears] ACt=OutPut[A+14,]';
  
  vector[nyears] JCtProp=JCt ./ (JCt+ACt);

  
  vector[nyears] SSBStatus=OutPut[A+15,]';
  
  vector[nyears] availBStatus=OutPut[A+16,]';
  
  
  
  //trajectories under 20% or over 40%
  
  vector[nyears] TB40=rep_vector(0, nyears);
  vector[nyears] TB20=rep_vector(0, nyears);
  
  vector[nyears] AB40=rep_vector(0, nyears);
  vector[nyears] AB20=rep_vector(0, nyears);
  
    
  vector[nyears] JB40=rep_vector(0, nyears);
  vector[nyears] JB20=rep_vector(0, nyears);


  vector[nyears] SSB40=rep_vector(0, nyears);
  vector[nyears] SSB20=rep_vector(0, nyears);
  

  
  }

