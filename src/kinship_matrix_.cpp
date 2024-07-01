#include <Rcpp.h>
#include "Mozza.h"

//[[Rcpp::export]]
NumericMatrix kinship_matrix_(List L) {
  // il faut créer un vecteur de pointeurs... mozza::kinship_matrix ne fonctionne pas avec une liste.
  std::vector<mozza::zygote *> Z;
  for(auto a : L) {
    Rcpp::XPtr<mozza::zygote> pz = a;
    Z.push_back(pz);
  }
  return mozza::kinship_matrix(Z);
}

/*



Rcpp::XPtr<mozza::zygote> zygote(int ntiles, double mean_length_tiles = 20) {
  mozza::zygote * pz (new mozza::zygote(mozza::human_autosomes_b37, ntiles, mean_length_tiles));
  Rcpp::XPtr<mozza::zygote> xpz(pz, true);
  return xpz;
}

//[[Rcpp::export]]
Rcpp::XPtr<mozza::zygote> zygote_probs(NumericVector probaTiles, double mean_length_tiles = 20) {
  mozza::zygote * pz (new mozza::zygote(mozza::human_autosomes_b37, probaTiles, mean_length_tiles));
  Rcpp::XPtr<mozza::zygote> xpz(pz, true);
  return xpz;
}

*/
