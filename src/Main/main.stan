// --------------------------------------------------------------
// Length-Based Age-Structured Model by Kyuhan Kim
// Copyright © 2023 Kyuhan Kim. All rights reserved.
// Contact: kh2064@gmail.com for questions
// MIT License: https://opensource.org/licenses/MIT
// --------------------------------------------------------------

functions{
  #include "../Submodel/Selex.stan"
  #include "../Submodel/StockRecruitment.stan"
  #include "../Submodel/eqQuant.stan"
  #include "../Submodel/dynamics.stan"
  #include "../Submodel/HelperFunctions.stan"
  }

data {
  int nbins;        // number of length bins
  int nCatchYears;  // number of years of catch data
  int A;            // number of age classes
  int nProj;        // number of projection years
  int nyears;       // number of years
  int nplatoons;    // number of platoons
  int FindBrps;     // find Biological Reference Points?
  
  vector[nbins] bins;          // length bins
  vector[nCatchYears] JYt;     // juvenile catch data
  vector[nCatchYears] AYt;     // adult catch data
  vector[nbins] LengthWeight;  // length-weight relationship
  vector[nbins] Maturity;      // maturity-at-length
  vector[nyears] Release;      // release numbers
  vector[nplatoons] Rprop;     // proportion of recruits in each platoon
  
  real FemaleProp;             // female proportion
  
  real LrefForM;               // reference length for natural mortality
  
  array[nplatoons] matrix[A, nbins] LengthDist;  // length-at-age distribution
  
  real steepPoint;             // steepness point (default: 0.2)
  //real gamma;
  
  real NatExponent;            // length-dependent natural mortality model exponent
  
  //real L05;
  //real L95;
  real sigR;                   // recruitment standard deviation
  real JuvMaxVul;              // maximum juvenile vulnerability
  real AduMaxVul;              // maximum adult vulnerability
              
  real L05_low;                // lower bound of L05
  real L05_up;                 // upper bound of L05
  
  real L95_low;               // lower bound of L95
  real L95_up;                // upper bound of L95
  
  }

parameters {
  
  real R0;                     // unfished recruitment
  real Natural_M;              // natural mortality
  real steepness;              // steepness
  real gamma;                  // depensation parameter
  real L95;                    // length at 95% selectivity
  real L05;                    // length at 5% selectivity
  vector[nyears+A-1] Rdev;     // recruitment deviations
  vector[nyears] L05dev;       // L05 deviations
  vector[nyears] L95dev;       // L95 deviations
  
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

