//[[Rcpp::depends(gaston)]]
#include <Rcpp.h>
#include "Mozza.h"
#include <fstream>

using namespace Rcpp;

// Z = une liste (pointeurs vers des zygotes)
//[[Rcpp::export]]
void drop_genotypes_to_vcf(List Z, std::string filename, XPtr<matrix4> Haplos, IntegerVector chr, NumericVector dist, CharacterVector id, 
                           IntegerVector pos, CharacterVector A1, CharacterVector A2) {

  // open file for append
  std::ofstream output(filename, std::ios_base::app); 
  std::vector<mozza::zygote *> ZYG;
  for(auto a : Z) {
    Rcpp::XPtr<mozza::zygote> pz = a;
    ZYG.push_back(pz);
  }

  mozza::mappedBed<IntegerVector, NumericVector, CharacterVector> MB(Haplos, chr, dist, id, pos, A1, A2);
  zygote_to_vcf(ZYG, MB, output);
}
