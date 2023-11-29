//Equilibrium N
matrix EqN(vector Nat_i,
           int nbins,
           int nplatoons,
           real R0,
           array[] matrix LengthDist,
           vector Rprop,
           int A) {
             
  matrix[nplatoons, A] unfishedN=rep_matrix(0., nplatoons, A);

  
  for (a in 1:A) {
      
      if(a==1) {
        for(p in 1:nplatoons) {
        unfishedN[p,a]=R0*Rprop[p];
        }
      } else if (a>1 || a<A) {
        for (p in 1:nplatoons) {
        for (i in 1:nbins) {
        unfishedN[p,a]+=unfishedN[p,a-1]*LengthDist[p,a-1, i]*exp(-Nat_i[i]);
        }
        }
      } else if (a==A) {
        for (p in 1:nplatoons) {
        for (i in 1:nbins) {
          unfishedN[p,a]+=unfishedN[p,a-1]*LengthDist[p,a-1, i]*(exp(-Nat_i[i])/(1.-exp(-Nat_i[i])));
        }
        }
      }
      
      }
    
    return unfishedN;
    
    }
    