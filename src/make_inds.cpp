//[[Rcpp::depends(gaston)]]
#include <Rcpp.h>
#include "mozza.h"
#include "gaston/matrix4.h"
using namespace Rcpp;

// cette fonction fait des individus indépendants
// avec des haplos mosaiques avec tuiles de longueurs 
// length_tiles parmi n_haps haplotypes
// Ils sont renvoyés comme un vecteur de mozza::zygote
std::vector<mozza::zygote> make_inds(int n, int n_haps, double length_tiles = 20.) {
  std::vector<mozza::zygote> ZYG;
  for(int i = 0; i < n; i++) {
    // mozza::mosaic H1(mozza::human_autosomes_b37, n_haps, length_tiles);
    // mozza::mosaic H2(mozza::human_autosomes_b37, n_haps, length_tiles);
    // mozza::zygote Z( H1, H2 );
    ZYG.push_back( mozza::zygote(mozza::human_autosomes_b37, n_haps, length_tiles) );
  }
  return ZYG;
}

//[[Rcpp::export]]
List make_inds(int n, double length_tiles, XPtr<matrix4> Haplos, IntegerVector chr, NumericVector dist, bool kinship = false) {
  int n_haps = Haplos->ncol; // chaque haplotype = un "individu"
  std::vector<mozza::zygote> ZYG { make_inds(n, n_haps, length_tiles) }; 
  
  List L;
  L["bed"] = drop_to_bed_matrix(ZYG, Haplos, chr, dist);
  if(kinship) 
    L["kinship"] = kinship_matrix(ZYG);
  return L;
}
