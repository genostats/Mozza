//[[Rcpp::depends(gaston)]]
#include <Rcpp.h>
#include "mozza.h"
using namespace Rcpp;

//[[Rcpp::export]]
List cousins_1stdegree(int n, double tile_length, XPtr<matrix4> Haplos, IntegerVector chr, NumericVector dist, 
                       bool kinship = false, bool fraternity = false) {
  std::vector<mozza::zygote> ZYG; 
  int n_haps = Haplos->ncol; // chaque haplotype = un "individu"
  
  for(int i = 0; i < n; i++) {
    mozza::zygote GM(mozza::human_autosomes_b37, n_haps, tile_length);
    mozza::zygote GF(mozza::human_autosomes_b37, n_haps, tile_length);
    mozza::zygote A = GM + GF;
    mozza::zygote B = GM + GF;
    mozza::zygote A1(mozza::human_autosomes_b37, n_haps, tile_length);
    mozza::zygote B1(mozza::human_autosomes_b37, n_haps, tile_length);
    ZYG.push_back(A + A1);
    ZYG.push_back(B + B1);
  }
  List L;
  L["bed"] = drop_to_bed_matrix(ZYG, Haplos, chr, dist);
  if(kinship) 
    L["kinship"] = kinship_matrix(ZYG);
  if(fraternity) 
    L["fraternity"] = fraternity_matrix(ZYG);
  return L;
}
