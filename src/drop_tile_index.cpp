#include <Rcpp.h>
#include "Mozza.h"
#include "Mozza/drop_tile_index.h"

//[[Rcpp::export]]
Rcpp::IntegerVector haplo_drop_tile_index(Rcpp::XPtr<mozza::mosaic> xpm, IntegerVector chr, NumericVector pos) {
  std::vector<int> tile_index;
  tile_index.reserve(chr.size());
  drop_tile_index(*xpm, chr, pos, true, tile_index);
  return wrap(tile_index);
}

//[[Rcpp::export]]
Rcpp::List zygote_drop_tile_index(Rcpp::XPtr<mozza::zygote> xpz, IntegerVector chr, NumericVector pos) {
  std::vector<int> tile_index1;
  tile_index1.reserve(chr.size());
  drop_tile_index(xpz->first, chr, pos, true, tile_index1);

  std::vector<int> tile_index2;
  tile_index2.reserve(chr.size());
  drop_tile_index(xpz->second, chr, pos, true, tile_index2);

  Rcpp::List L(2);
  L[0] = wrap(tile_index1);
  L[1] = wrap(tile_index2);
  return L;
}


