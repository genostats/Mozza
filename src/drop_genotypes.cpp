//[[Rcpp::depends(gaston)]]
#include <Rcpp.h>
#include "Mozza.h"
using namespace Rcpp;

// Z = une liste (pointeurs vers des zygotes)
//[[Rcpp::export]]
List drop_genotypes(List Z, XPtr<matrix4> Haplos, IntegerVector chr, NumericVector dist, bool phased) {
  std::vector<mozza::zygote *> ZYG;
  for(auto a : Z) {
    Rcpp::XPtr<mozza::zygote> pz = a;
    ZYG.push_back(pz);
  }

  mozza::mappedBed<IntegerVector, NumericVector> MB(Haplos, chr, dist);
  List L;
  L["bed"] = zygote_to_bed_matrix(ZYG, MB, phased);
  return L;
}
