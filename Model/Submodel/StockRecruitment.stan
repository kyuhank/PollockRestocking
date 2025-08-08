//Stock-recruitment relationship
real StockRecruitment(real R0,
                      real steepness,
                      real SSB,
                      real SSB0,
                      real gamma,
                      real steepPoint) {
                        
      real alpha;
      real beta;
      real Rec;
                      
      alpha=(steepness*R0*(1. -steepPoint^(gamma)))/(SSB0^(gamma)*steepPoint^(gamma)*(1.-steepness));
      beta=(1./(steepPoint^(gamma))*steepness-1.)/(SSB0^(gamma)*(1. -steepness));                  
      Rec=(alpha*SSB^(gamma))/(1.+beta*SSB^(gamma));
      
      return Rec;                  
      
      }