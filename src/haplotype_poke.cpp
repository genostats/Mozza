#include <Rcpp.h>
#include "Mozza.h"

// remplace la tuile dont l'indice (en partant du premier chromosome) est donné dans tile_index
// par la valeur données dans tile
// RIndex = true si les index commencent à 1
// la fonction suppose que tile_index est trié dans l'ordre croissant
template<typename intVec>
void haplotype_poke(mozza::mosaic & M, intVec tile_index, intVec tile, bool RIndex = true) {
  if(M.chrs == 0) throw std::runtime_error("Empty haplotype");
  int n = tile_index.size();
  if(tile.size() != n) throw std::runtime_error("tile and tile_index should have the same size");
  int chr = 0;
  unsigned int offset = 0;
  unsigned int ntiles = M.tiles[0].size();
  const int ri = RIndex?1:0;
  for(int a = 0; a < n; a++) {
    unsigned int  i = tile_index[a] - ri;
    unsigned int ti = tile[a];
    while(i >= offset + ntiles) {
      offset += ntiles;
      if(++chr > M.chrs) throw std::runtime_error("Bad tile index");
      ntiles = M.tiles[chr].size();
    }
    // so (i - offset) < ntiles
    M.tiles[chr][i - offset] = ti;
  }
}



//[[Rcpp::export]]
void haplotype_poke(Rcpp::XPtr<mozza::mosaic> xpm, IntegerVector tile_index, IntegerVector tile) {
  haplotype_poke(*xpm, tile_index, tile);
}


//[[Rcpp::export]]
void zygote_poke(Rcpp::XPtr<mozza::zygote> xpz, IntegerVector tile_index1, IntegerVector tile1, IntegerVector tile_index2, IntegerVector tile2) {
  haplotype_poke(xpz->first, tile_index1, tile1);
  haplotype_poke(xpz->second, tile_index2, tile2);
}

