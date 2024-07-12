#include <Rcpp.h>
#include "Mozza.h"

//[[Rcpp::export]]
Rcpp::XPtr<mozza::mosaic> haplotype(int ntiles, double mean_length_tiles = 20) {
  mozza::mosaic * pm (new mozza::mosaic(mozza::Autosomes(), ntiles, mean_length_tiles));
  Rcpp::XPtr<mozza::mosaic> xpm(pm, true);
  return xpm;
}

//[[Rcpp::export]]
Rcpp::XPtr<mozza::mosaic> haplotype_probs(NumericVector probaTiles, double mean_length_tiles = 20) {
  mozza::mosaic * pm (new mozza::mosaic(mozza::Autosomes(), probaTiles, mean_length_tiles));
  Rcpp::XPtr<mozza::mosaic> xpm(pm, true);
  return xpm;
}

