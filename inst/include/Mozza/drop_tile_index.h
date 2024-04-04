#include "Mozza/mosaic.h"
#include "Mozza/zygote.h"

#ifndef _mozza_drop_tile_index_
#define _mozza_drop_tile_index_

namespace mozza {
// pushes back in tile_index the tile index (numbering from first chromosome) of variants with position
// given in chr and pos (sorted in ascending order
// RIndex = true : starts numbering at 1 (applies to chromosomes and tile indices !)
template<typename intVec, typename floatVec>
void drop_tile_index(mosaic & M, intVec chr, floatVec pos, bool RIndex, std::vector<int> & tile_index) {
  int n = chr.size();
  if(pos.size() != n) throw std::runtime_error("chr and pos should have the same size");
  const int ri = RIndex?1:0;
  M.set_cursor(); // position cursor at beginning of chr 0
  unsigned int current_chr = 0;
  unsigned int offset = 0;
  unsigned int ntiles = M.tiles[0].size();
  for(int a = 0; a < n; a++) {
    unsigned int c = chr[a] - ri;
    auto p = pos[a];
    // change chromosome if needed
    while(current_chr < c) {
      offset += ntiles;
      current_chr++;
      M.set_cursor(current_chr);    
      ntiles = M.tiles[current_chr].size();
    }
    // current_chr == c
    M.set_cursor(current_chr, p);
    tile_index.push_back(offset + M.i_cursor + ri);
  }
}

} // end namespace

#endif
