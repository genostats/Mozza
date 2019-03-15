#include <Rcpp.h>
#include "mosaic.h"
using namespace Rcpp;

namespace mozza {

void mosaic::print_chr(int i) {
  Rcout << "  chr " <<  i << ": " << tiles[i].size() << " tile(s) ";
  for(int j = 0; j < tiles[i].size() ; j++) {
    Rcout << tiles[i][j] << " [" << bpoints[i][j] << "] ";
  }
  Rcout << "\n";
}

void mosaic::print() {
  for(int i = 0; i < chrs; i++) {
    print_chr(i);
  }
}

void mosaic::print_cursor() {
  Rcout << "cursor at " << cursor_chr << ": " << cursor_pos;
  Rcout << " (i = " << i_cursor << " tile = " << tile_at_cursor() << ")\n";
}
};
