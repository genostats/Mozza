#include <Rcpp.h>
#include "Mozza.h"
#include "Mozza/drop_tile_index.h"


Rcpp::List group_by_tile_index(std::vector<int> & tile_index) {
  Rcpp::List L;
  int n = tile_index.size();
  if(n == 0) return L;

  std::vector<int> group;
  int tile = tile_index[0];
  for(int i = 0; i < n; i++) {
    if(tile_index[i] != tile) {
      L[ std::to_string(tile) ] = wrap(group);
      group.clear();
      tile = tile_index[i];
    }
    group.push_back(i + 1);
  }
  L[ std::to_string(tile) ] = wrap(group);
  return L;
}

//[[Rcpp::export]]
Rcpp::List haplo_group_by_tile_index(Rcpp::XPtr<mozza::mosaic> xpm, IntegerVector chr, NumericVector pos) {
  std::vector<int> tile_index;
  tile_index.reserve(chr.size());
  drop_tile_index(*xpm, chr, pos, true, tile_index);
  return group_by_tile_index(tile_index);
}


//[[Rcpp::export]]
Rcpp::List zygote_group_by_tile_index(Rcpp::XPtr<mozza::zygote> xpz, IntegerVector chr, NumericVector pos) {
  std::vector<int> tile_index1;
  tile_index1.reserve(chr.size()); 
  drop_tile_index(xpz->first, chr, pos, true, tile_index1);

  std::vector<int> tile_index2;
  tile_index2.reserve(chr.size());
  drop_tile_index(xpz->second, chr, pos, true, tile_index2);

  Rcpp::List L(2);
  L[0] = group_by_tile_index(tile_index1);
  L[1] = group_by_tile_index(tile_index2);
  return L;
}

