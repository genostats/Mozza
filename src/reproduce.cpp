
#include <Rcpp.h>
#include "Mozza.h"

//[[Rcpp::export]]
Rcpp::XPtr<mozza::zygote> reproduce(Rcpp::XPtr<mozza::zygote> z1, Rcpp::XPtr<mozza::zygote> z2) {
  mozza::zygote * pzo (new mozza::zygote(*z1, *z2));
  Rcpp::XPtr<mozza::zygote> xpzo(pzo, true);
  return xpzo;
}
