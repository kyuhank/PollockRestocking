// --------------------------------------------------------------
// Length-Based Age-Structured Model by Kyuhan Kim
// Copyright © 2023 Kyuhan Kim. All rights reserved.
// Contact: kh2064@gmail.com for questions
// MIT License: https://opensource.org/licenses/MIT
// --------------------------------------------------------------

//dynamics loop
matrix DynamicsN(int nCatchYears,
                 int A,
                 int nplatoons,
                 vector bins,
                 vector LengthWeight,
                 vector Nat_i,
                 vector Maturity,
                 vector JVulMean,
                 matrix JVul,
                 vector AVul,
                 real R0,
                 real sigR,
                 real steepness,
                 real gamma,
                 real FemaleProp,
                 array[] matrix LengthDist,
                 vector JYt,
                 vector AYt,
                 real steepPoint,
                 int nbins,
                 vector Release,
                 vector Rdev,
                 vector Rprop,
                 int nProj,
                 int CatchSepYr) {
              
  
  int nyears=nCatchYears+nProj;

  array[nplatoons] matrix[A, nyears] Npat=rep_array(rep_matrix(0., A, nyears), nplatoons);
  array[nplatoons] matrix[A, A] preN=rep_array(rep_matrix(0., A, A), nplatoons);
  
  matrix[A, nyears] JCat=rep_matrix(0., A, nyears);
  matrix[A, nyears] ACat=rep_matrix(0., A, nyears);
  
  vector[nyears] Bt=rep_vector(0., nyears);
  vector[nyears] availBt=rep_vector(0., nyears);
  vector[nyears] SSBt=rep_vector(0., nyears);
  vector[nyears] TJBt=rep_vector(0., nyears);
  vector[nyears] TABt=rep_vector(0., nyears);
  
  vector[nyears] JBt=rep_vector(0., nyears); 
  vector[nyears] ABt=rep_vector(0., nyears);
  
  
  vector[nyears] SSBrel=rep_vector(0., nyears);
  vector[nyears] JBrel=rep_vector(0., nyears); 
  vector[nyears] ABrel=rep_vector(0., nyears);
  
  
  
  vector[nyears] JCt=rep_vector(0., nyears);
  vector[nyears] ACt=rep_vector(0., nyears);
  
  vector[nyears] JHt;
  vector[nyears] AHt;
  
  vector[nyears] BStatus;
  vector[nyears] JBStatus;
  vector[nyears] ABStatus;
  vector[nyears] SSBStatus;
  vector[nyears] availBStatus;
  
  vector[nyears] JuvProp;
  
  matrix[A+19, nyears] OutPut; 
  
  real SSB0=0.;
  real B0=0.;
  real availB0=0.;
  real JB0=0.;
  real AB0=0.;
  real SSBref=0.;
  real ABref=0.;
  real JBref=0.;
  real TJB0=0.;
  real TAB0=0.;
  
  vector[A] preSSBt=rep_vector(0., A);
  
  matrix[nplatoons,A] unfishedN=EqN(Nat_i,
                                    nbins,
                                    nplatoons,
                                    R0,
                                    LengthDist,
                                    Rprop,
                                    A);  
    
    
 //Equilibrium quantities    
 for(a in 1:A) {
   for (p in 1:nplatoons) {
    for(i in 1:nbins) {
    SSB0+=unfishedN[p,a]*LengthDist[p,a,i]*Maturity[i]*FemaleProp*LengthWeight[i];
    B0+=unfishedN[p,a]*LengthDist[p,a,i]*LengthWeight[i];
    JB0+=unfishedN[p,a]*LengthDist[p,a,i]*LengthWeight[i]*JVulMean[i];
    AB0+=unfishedN[p,a]*LengthDist[p,a,i]*LengthWeight[i]*AVul[i];
    TJB0+=unfishedN[p,a]*LengthDist[p,a,i]*LengthWeight[i]*(1.-Maturity[i]);
    TAB0+=unfishedN[p,a]*LengthDist[p,a,i]*LengthWeight[i]*Maturity[i];
    }
  }
 }
    availB0=JB0+AB0;
        
      
    
    for (p in 1:nplatoons) {  
      preN[p,1:A,1]=to_vector(unfishedN[p,1:A]);
    }              
                
 //for initial N with random Rdev                        
 for(j in 2:(A)) {
   
   for (p in 1:nplatoons) {
     for (a in 1:A) {
       for (i in 1:nbins) {
       preSSBt[j-1]+=preN[p,a,j-1]*LengthDist[p,a,i]*Maturity[i]*LengthWeight[i]*FemaleProp;
       }
     }
   }
   
  for (a in 1:A) {
    
    if(a==1) {
      for (p in 1:nplatoons) {
        preN[p,a,j]=Rprop[p]*exp(log(StockRecruitment(R0, steepness, preSSBt[j-1], SSB0, gamma, steepPoint))+Rdev[j-1]-0.5*sigR^2);
      }
    } else if (a>1 || a<A) {
      for (p in 1:nplatoons) {
        for (i in 1:nbins) {
          preN[p,a,j]+=preN[p,a-1,j-1]*LengthDist[p,a-1, i]*exp(-Nat_i[i]);
        }
      } 
    } else if (a==A) {
      for (p in 1:nplatoons) {
        for (i in 1:nbins) {
          preN[p,a,j]+=preN[p,a-1,j-1]*LengthDist[p,a-1, i]*exp(-Nat_i[i])+preN[p,a,j-1]*LengthDist[p,a, i]*exp(-Nat_i[i]);
        }
      }
    }
    
  }
  
 }
  

  Npat[1:nplatoons,1:A,1]=preN[1:nplatoons,1:A, A];

                
  for(a in 1:A) {
    for (p in 1:nplatoons) {
    for (i in 1:nbins) {
      Bt[1]+=Npat[p,a,1]*LengthDist[p,a,i]*LengthWeight[i];
      JBt[1]+=Npat[p,a,1]*LengthDist[p,a,i]*LengthWeight[i]*JVul[i,1];
      ABt[1]+=Npat[p,a,1]*LengthDist[p,a,i]*LengthWeight[i]*AVul[i];
      SSBt[1]+=Npat[p,a,1]*LengthDist[p,a,i]*Maturity[i]*LengthWeight[i]*FemaleProp;
      TJBt[1]+=Npat[p,a,1]*LengthDist[p,a,i]*LengthWeight[i]*(1.-Maturity[i]);
      TABt[1]+=Npat[p,a,1]*LengthDist[p,a,i]*LengthWeight[i]*Maturity[i];
    }
  }
  }
  
  availBt[1]=JBt[1]+ABt[1];
  
  BStatus[1]=Bt[1]/B0;
  JBStatus[1]=JBt[1]/JB0;
  ABStatus[1]=ABt[1]/AB0;
  SSBStatus[1]=SSBt[1]/SSB0;
  availBStatus[1]=availBt[1]/availB0;
  JuvProp[1]=(JBt[1]/(JBt[1]+ABt[1]));

  for(t in 2:nyears) {
    
    if (t <= (nCatchYears+1) ) {
    JHt[t-1]=JYt[t-1]/JBt[t-1];
    AHt[t-1]=AYt[t-1]/ABt[t-1];
    } else  {
    JHt[t-1]=0.;
    AHt[t-1]=0.;
    }
    
    if( CatchSepYr>0 && t<=(CatchSepYr+1)) {
      JHt[t-1]=((JYt[t-1]+AYt[t-1])*JuvProp[t-1])/JBt[t-1];
      AHt[t-1]=((JYt[t-1]+AYt[t-1])*(1-JuvProp[t-1]))/ABt[t-1];
    }
    
  
  for(a in 1:A){
    
    if(a==1){
      
      for (p in 1:nplatoons) {
      Npat[p,a,t]=exp(log(StockRecruitment(R0, steepness, SSBt[t-1], SSB0, gamma, steepPoint))+Rdev[A+(t-1)]-0.5*sigR^2)*Rprop[p];
      if(p==nplatoons) {
      Npat[p,a,t]=Npat[p,a,t]+Release[t];
      }
      }

    } else if (a>1 || a<A) {
      for (p in 1:nplatoons) {
     for(i in 1:nbins) {
         Npat[p,a,t]+=exp(log( Npat[p,a-1, t-1]*LengthDist[p,a-1,i]*(1.-(AHt[t-1]*AVul[i]+JHt[t-1]*JVul[i,t-1]))*exp(-Nat_i[i]) ));
     }
      }
    } else if (a==A) {
      for (p in 1:nplatoons) {
      for(i in 1:nbins) {
        Npat[p,a,t]+=exp(log(  Npat[p,a-1, t-1]*LengthDist[p,a-1,i]*(1.-(AHt[t-1]*AVul[i]+JHt[t-1]*JVul[i,t-1]))*exp(-Nat_i[i]) +  Npat[p,a, t-1]*LengthDist[p,a,i]*(1.-(AHt[t-1]*AVul[i]+JHt[t-1]*JVul[i,t-1]))*exp(-Nat_i[i]) ));
      }
    }
    }
         
  }
  
  for(a in 1:A) {
    for (p in 1:nplatoons) {
    for(i in 1:nbins) {
    SSBt[t]+=Npat[p,a,t]*LengthDist[p,a,i]*Maturity[i]*LengthWeight[i]*FemaleProp;
    Bt[t]+=Npat[p,a,t]*LengthDist[p,a,i]*LengthWeight[i];
    
    JBt[t]+=Npat[p,a,t]*LengthDist[p,a,i]*LengthWeight[i]*JVul[i,t];
    ABt[t]+=Npat[p,a,t]*LengthDist[p,a,i]*LengthWeight[i]*AVul[i];
    
    TJBt[t]+=Npat[p,a,t]*LengthDist[p,a,i]*LengthWeight[i]*(1-Maturity[i]);
    TABt[t]+=Npat[p,a,t]*LengthDist[p,a,i]*LengthWeight[i]*(Maturity[i]);
    
    
   }
    }
  }
  
  availBt[t]=ABt[t]+JBt[t];
  
  BStatus[t]=Bt[t]/B0;
  JBStatus[t]=TJBt[t]/TJB0;
  ABStatus[t]=TABt[t]/TAB0;
  SSBStatus[t]=SSBt[t]/SSB0;
  availBStatus[t]=availBt[t]/availB0;
  JuvProp[t]=(JBt[t]/(JBt[t]+ABt[t]));

  }
  
  for(t in 1:nyears){
   for(a in 1:A){
     for (p in 1:nplatoons) {
     for(i in 1:nbins) {
      JCat[a,t]+=Npat[p,a, t]*LengthDist[p,a,i]*(JHt[t]*JVul[i,t]);
      ACat[a,t]+=Npat[p,a, t]*LengthDist[p,a,i]*(AHt[t]*AVul[i]);
      }
     }
    JCt[t]+=JCat[a,t];
    ACt[t]+=ACat[a,t];
    }
  }
  
  
  
  matrix[A, nyears] Nat=rep_matrix(0., A, nyears);
  
  for(t in 1:nyears){
    for (p in 1:nplatoons) {
      for(a in 1:A){
        Nat[a,t]+=Npat[p,a,t];
      }
    }
  }
  
  SSBref=mean(SSBt[25:35]);
  ABref=mean(TABt[25:35]);
  JBref=mean(TJBt[25:35]);
  
  SSBrel=SSBt/SSBref;
  ABrel=TABt/ABref;
  JBrel=TJBt/JBref;
  
  
  OutPut[1:A,]=Nat;
  OutPut[A+1,]=TJBt';
  OutPut[A+2,]=TABt';
  OutPut[A+3,]=SSBt';
  OutPut[A+4,]=Bt';
  OutPut[A+5,]=BStatus';
  OutPut[A+6,]=JBStatus';
  OutPut[A+7,]=ABStatus';
  OutPut[A+8,]=JHt';
  OutPut[A+9,]=AHt';
  OutPut[A+10,]=availBt';
  OutPut[A+11,1]=SSB0;
  OutPut[A+11,2]=availB0;
  OutPut[A+11,3]=B0;
  OutPut[A+11,4]=JB0;
  OutPut[A+11,5]=AB0;
  OutPut[A+12,]=(Release ./ (Nat[1,1:nyears]'+Release))';
  OutPut[A+13,]=JCt';
  OutPut[A+14,]=ACt';
  OutPut[A+15,]=SSBStatus';
  OutPut[A+16,]=availBStatus';
  OutPut[A+17,]=SSBrel';
  OutPut[A+18,]=JBrel';
  OutPut[A+19,]=ABrel';

  return OutPut;              
                  
  }