
#include <Rcpp.h>
#include "Mozza.h"

//[[Rcpp::export]]
Rcpp::XPtr<mozza::zygote> reproduce_(Rcpp::XPtr<mozza::zygote> z1, Rcpp::XPtr<mozza::zygote> z2) {
  mozza::zygote * pzo (new mozza::zygote(*z1, *z2));
  Rcpp::XPtr<mozza::zygote> xpzo(pzo, true);
  return xpzo;
}

//[[Rcpp::export]]
Rcpp::List reproduce_vec(Rcpp::List Z1, Rcpp::List Z2) {
  unsigned int n = Z1.size();
  if(n != Z2.size()) stop("Not the same length");
  Rcpp::List ZYG(n);
  for(unsigned int i = 0; i < n; i++) {
    ZYG[i] = reproduce_(Z1[i], Z2[i]);
  }
  ZYG.attr("class") = CharacterVector::create("zygote");
  return ZYG;
}
