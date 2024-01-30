#include <Rcpp.h>
#include "mosaic.h"

namespace mozza {

// longueur IBD entre haplotypes M1 et M2
// note : on ne vérifie pas que la longueur des chromosomes / la fin de la dernière tuile coincident...
// NOTE la position du curseur est modifiée (on pourrait la sauver et remettre le curseur à sa place ?)
double IBD_sharing(mosaic & M1, mosaic & M2) {
  double shared = 0.;
  for(int i = 0; i < M1.chrs; i++)  {
    std::vector<double> pos(M1.bpoints[i].size()+M2.bpoints[i].size());
    std::merge(M1.bpoints[i].begin(), M1.bpoints[i].end(), 
               M2.bpoints[i].begin(), M2.bpoints[i].end(),
               pos.begin());
    double a = 0.;
    M1.set_cursor(i);
    M2.set_cursor(i);
    for(auto & b : pos) {
      M1.forward_cursor(b);
      M2.forward_cursor(b);
      if(M1.tile_at_cursor() == M2.tile_at_cursor())
        shared += (b-a);
      a = b;
    }
  }
  return shared;
}

}
