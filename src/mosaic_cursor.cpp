#include <Rcpp.h>
#include "mosaic.h"
using namespace Rcpp;

namespace mozza {
 
void mosaic::set_cursor(unsigned int chr, double pos){
  if(chr < chrs) {
    cursor_chr = chr;
    if(0 <= pos && pos <= chr_len[chr]) {
      cursor_pos = pos;
    } else {
      Rcpp::warning("Trying to set cursor pos outside of chr");
      cursor_pos = chr_len[chr];
    }
  } else
    stop("Trying to set cursor to an inexistent chr");
  i_cursor = 0;
  // on incrémente i_cursor tant que le bpoint de la tile n'a pas dépassé pos
  while(bpoints[chr][i_cursor] < pos)
    i_cursor++;
}


void mosaic::step_cursor(double len) {
  cursor_pos += len;
  if(cursor_pos > chr_len[cursor_chr]) 
    cursor_pos = chr_len[cursor_chr];
  // on incrémente i_cursor tant que le bpoint de la tile n'a pas dépassé pos
  while(bpoints[cursor_chr][i_cursor] < cursor_pos)
    i_cursor++;
}


void mosaic::forward_cursor(double pos) {
  if(pos >= cursor_pos) // on avance si c'est plus loin sinon on reste sur place
    cursor_pos = pos;
  else {
    Rcpp::warning("Trying to forward cursor to a 'backward' position " + std::to_string(cursor_chr) + ":" + std::to_string(cursor_pos) + " to " + std::to_string(pos));
    return;
  }
  
  if(cursor_pos > chr_len[cursor_chr]) 
    cursor_pos = chr_len[cursor_chr];
  // on incrémente i_cursor tant que le bpoint de la tile n'a pas dépassé pos
  while(bpoints[cursor_chr][i_cursor] < cursor_pos)
    i_cursor++;
}

void mosaic::step_cursor(double len, std::vector<int> & til, std::vector<double> & bpt) {
  cursor_pos += len;
  if(cursor_pos > chr_len[cursor_chr])
    cursor_pos = chr_len[cursor_chr];
  // on incrémente i_cursor tant que le bpoint de la tile n'a pas dépassé pos
  while(bpoints[cursor_chr][i_cursor] < cursor_pos) {
    til.push_back(tiles[cursor_chr][i_cursor]);
    bpt.push_back(bpoints[cursor_chr][i_cursor]);
    i_cursor++;
  }
}

// numéro de la tuile au curseur
// (les tuiles sont fermées à droite)
// int mosaic::tile_at_cursor const {
//   return tiles[cursor_chr][i_cursor];
// }

}
