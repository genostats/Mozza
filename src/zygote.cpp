#include <Rcpp.h>
#include "Mozza.h"

//[[Rcpp::export]]
Rcpp::XPtr<mozza::zygote> zygote_(int ntiles, double mean_length_tiles = 20) {
  mozza::zygote * pz (new mozza::zygote(mozza::Autosomes(), ntiles, mean_length_tiles));
  Rcpp::XPtr<mozza::zygote> xpz(pz, true);
  return xpz;
}

//[[Rcpp::export]]
Rcpp::XPtr<mozza::zygote> zygote_probs(NumericVector probaTiles, double mean_length_tiles = 20) {
  mozza::zygote * pz (new mozza::zygote(mozza::Autosomes(), probaTiles, mean_length_tiles));
  Rcpp::XPtr<mozza::zygote> xpz(pz, true);
  return xpz;
}


