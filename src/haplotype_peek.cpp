#include <Rcpp.h>
#include "Mozza.h"

Rcpp::List haplotype_peek(mozza::mosaic & M) {
  std::vector<int> chr;
  std::vector<double> pos;
  std::vector<int> tile;

  for(int i = 0; i < M.chrs; i++) {
    int ntile = M.tiles[i].size();
    for(int j = 0; j < ntile; j++) {
      tile.push_back( M.tiles[i][j] );
      pos.push_back( M.bpoints[i][j] );
      chr.push_back(i);
    }
  }   
  List L;
  L["chr"] = chr;
  L["pos"] = pos;
  L["tile"] = tile;
  return L;   
}

 
//[[Rcpp::export]]
Rcpp::List haplotype_peek(Rcpp::XPtr<mozza::mosaic> xpm) {
  return haplotype_peek(*xpm);
}

//[[Rcpp::export]]
Rcpp::List zygote_peek(Rcpp::XPtr<mozza::zygote> xpz) {
  Rcpp::List L(2);
  L[0] = haplotype_peek(xpz->first);
  L[1] = haplotype_peek(xpz->second);
  return L;
}
