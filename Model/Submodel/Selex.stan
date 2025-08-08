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