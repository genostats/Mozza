//[[Rcpp::depends(gaston)]]
#include <Rcpp.h>
#include "mozza.h"
#include "gaston/matrix4.h"
using namespace Rcpp;

// cette fonction fait des paires d'individus avec des haplos 
// H1 H2, M1 M2, ou H1 et M1 sont indépendants, mais H2 et M2
// partagent des segments de longueur le1 entrecoupés par des
// segments de longueur le2 (en espérance)
// Ils sont pushé dans ZYGS
// On renvoie leur coefficient k = 2*phi
std::vector<double> make_pairs(int n, double le1, double le2, std::vector<mozza::zygote> & ZYG, int n_haps = 100, double length_tiles = 20.) {
  std::vector<double> kin;
  for(int i = 0; i < n; i++) {
    mozza::mosaic H1(mozza::human_autosomes_b37, n_haps, length_tiles);
    mozza::mosaic H2(mozza::human_autosomes_b37, n_haps, length_tiles);

    mozza::mosaic M1(mozza::human_autosomes_b37, n_haps, length_tiles);
    mozza::mosaic AA(mozza::human_autosomes_b37, n_haps, length_tiles);

    mozza::zygote Z1( H1, H2 );
    mozza::zygote Z2( M1, mozza::mosaic(H2, AA, le1, le2) );

    auto IBD = IBD_length( Z1, Z2 );
    kin.push_back( (0.5*std::get<1>(IBD) + std::get<2>(IBD)) / mozza::length_human_autosomes_b37  );
    ZYG.push_back(Z1);
    ZYG.push_back(Z2);
  }
  return kin;
}

//[[Rcpp::export]]
List make_pairs(int N, double le1, double le2, double length_tiles, XPtr<matrix4> Haplos, IntegerVector chr, NumericVector dist) {
  std::vector<mozza::zygote> x; 
  int n_haps = Haplos->ncol; // chaque haplotype = un "individu"
  List L;
  L["kin"] = make_pairs(N, le1, le2, x, n_haps, length_tiles);
  L["bed"] = drop_to_bed_matrix(x, Haplos, chr, dist);
  return L;
}


