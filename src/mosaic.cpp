#include <Rcpp.h>
#include "mosaic.h"
using namespace Rcpp;

namespace mozza {

// constructeur avec une tuile 'tile' constante
mosaic::mosaic(const std::vector<double> & chr_len, int tile) : 
chrs(chr_len.size()), chr_len(chr_len), genome_length(std::accumulate(chr_len.begin(), chr_len.end(), 0.)), 
tiles(chrs), bpoints(chrs), cursor_chr(0), cursor_pos(0.), i_cursor(0) {
  for(int i = 0; i < chrs; i++) {
    tiles[i].push_back(tile);
    bpoints[i].push_back(chr_len[i]);
  }
}

// constructeur avec des tuiles de 0 à ntiles-1, de longueur prise dans une exponentielle 
// avec longueur moyenne mean_length
mosaic::mosaic(const std::vector<double> & chr_len, int ntiles, double mean_length) : 
chrs(chr_len.size()), chr_len(chr_len), genome_length(std::accumulate(chr_len.begin(), chr_len.end(), 0.)), 
tiles(chrs), bpoints(chrs), cursor_chr(0), cursor_pos(0.), i_cursor(0) {
  for(int i = 0; i < chrs; i++) {
    double le = chr_len[i];
    double pos = 0; 
    while(true) {
      // tire et ajoute une tuile
      tiles[i].push_back(R::runif(0,1)*ntiles);
      // position de fin de cette tuile
      pos += R::rexp(mean_length);
      if(pos < le) {
        bpoints[i].push_back(pos);
      } else {
        bpoints[i].push_back(le);
        break;
      }
    }
  }
}  

// constructeur qui mélange M1 et M2, en prenant dans M1 des morceaux de longueur
// ~ exp(le1) et dans M2 des morceaux de longueur ~ exp(le2)
mosaic::mosaic(mosaic & M1, mosaic & M2, double le1, double le2) :
chrs(M1.chrs), chr_len(M1.chr_len), genome_length(M1.genome_length), 
tiles(chrs),  bpoints(chrs), cursor_chr(0), cursor_pos(0.), i_cursor(0) {
  double p1 = le1/(le1 + le2); // proba d'être sur M1.
  if(M2.chrs != chrs)
    stop("Two mosaic haplotypes with different chrs");
  for(int i = 0; i < chrs; i++) {
    M1.set_cursor(i);
    M2.set_cursor(i);
    
    double le = chr_len[i];
    if(le != M2.chr_len[i])
      stop("Chromosomes of different lengths");
    // on tire sur quel haplo on commence
    bool on_m1 = (R::runif(0,1) < p1);
    // on avance...
    while(M1.cursor_pos < chr_len[i]) {
      double le = on_m1?R::rexp(le1):R::rexp(le2);
      if(on_m1) {
        M1.step_cursor(le, tiles[i], bpoints[i]);
        M2.step_cursor(le);
      } else {
        M1.step_cursor(le);
        M2.step_cursor(le, tiles[i], bpoints[i]);
      }
      // il faut ajouter le point de crossover
      if(on_m1)
        tiles[i].push_back(M1.tile_at_cursor());
      else
        tiles[i].push_back(M2.tile_at_cursor());
      bpoints[i].push_back(M1.cursor_pos);
      
      on_m1 = !on_m1; // on change
    }
  }
}

}; // end namespace