#include <Rcpp.h>
#include "Mozza.h"

namespace mozza {

  static std::vector<double> autosomes = human_autosomes_b37;
  static double length_autosomes = length_human_autosomes_b37;

  std::vector<double> & Autosomes() {
    return autosomes;
  }  

  double lengthAutosomes() {
    return length_autosomes;
  }
}

//[[Rcpp::export]]
NumericVector getAutosomes() {
  return wrap(mozza::Autosomes());
}

//[[Rcpp::export]]
void setAutosomes(NumericVector A) {
  mozza::autosomes.clear();
  mozza::length_autosomes = 0;
  for(auto le : A) {
    mozza::autosomes.push_back(le);
    mozza::length_autosomes += le;
  }
}
