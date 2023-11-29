// --------------------------------------------------------------
// Length-Based Age-Structured Model by Kyuhan Kim
// Copyright © 2023 Kyuhan Kim. All rights reserved.
// Contact: kh2064@gmail.com for questions
// MIT License: https://opensource.org/licenses/MIT
// --------------------------------------------------------------



real bounded_inv_logit(real y, real Lower, real Upper) {
  
  real x=(Upper-Lower)/(1.0+exp(-y))+Lower;
  
  return x;
  
}


real bounded_logit(real x, real Lower, real Upper) {
  
  real y=-log(((Upper-Lower)/(x-Lower))-1);
  
  return y;
  
}




