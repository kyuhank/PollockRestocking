// --------------------------------------------------------------
// Length-Based Age-Structured Model by Kyuhan Kim
// Copyright © 2023 Kyuhan Kim. All rights reserved.
// Contact: kh2064@gmail.com for questions
// MIT License: https://opensource.org/licenses/MIT
// --------------------------------------------------------------


//Selectivity function    
vector Sel(vector bins,
           real L05, 
           real L95,
           int nbins
           ) {
    
    vector[nbins] Selec;
    
    for(i in 1:nbins) {
      Selec[i]=19./(19.+exp(-log(361.)*( (bins[i]-L95)/(L95-L05))));
    }
             
    return Selec;         
  }