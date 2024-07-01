
#include <Rcpp.h>
#include "Mozza.h"

//[[Rcpp::export]]
Rcpp::XPtr<mozza::zygote> reproduce(Rcpp::XPtr<mozza::zygote> xpz1, Rcpp::XPtr<mozza::zygote> xpz2) {
  mozza::zygote * pzo (new mozza::zygote(*xpz1, *xpz2));
  Rcpp::XPtr<mozza::zygote> xpzo(pzo, true);
  return xpzo;
}
