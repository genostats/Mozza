//[[Rcpp::depends(gaston)]]
#include <Rcpp.h>
#include "mozza.h"
using namespace Rcpp;

//[[Rcpp::export]]
List nuclear_families(int nb_fams, int nb_offsprings, double tile_length, XPtr<matrix4> Haplos, IntegerVector chr, NumericVector dist, bool kinship = false) {
  std::vector<mozza::zygote> ZYG; 
  int n_haps = Haplos->ncol; // chaque haplotype = un "individu"
  
  for(int i = 0; i < nb_fams; i++) {
    mozza::zygote M(mozza::human_autosomes_b37, n_haps, tile_length);
    mozza::zygote F(mozza::human_autosomes_b37, n_haps, tile_length);
    ZYG.push_back(M);
    ZYG.push_back(F);
    for(int j = 0; j < nb_offsprings; j++)
      ZYG.push_back(M+F);
  }
  List L;
  L["bed"] = drop_to_bed_matrix(ZYG, Haplos, chr, dist);
  if(kinship) 
    L["kinship"] = kinship_matrix(ZYG);
  return L;
}