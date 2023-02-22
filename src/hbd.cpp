#include <Rcpp.h>
#include "mosaic.h"
#include "zygote.h"
#include "HBD_at_point.h"
/* 
 *
 * calcul des longueurs HBD
 * écrit à partir des fonctions IBD_at_point et IBD_length (de ibd.cpp)
 *
 */


namespace mozza {

// longueur HBD des deux haplotypes du zygote
// NOTE la position des curseurs est modifiée (on pourrait la sauver et remettre le curseur à sa place ?)
//      cf ibd.cpp pour plus de réflexions sur la façon de rendre ça propre
double HBD_length(zygote & Z) {
  double HBD = 0;
  for(int i = 0; i < Z.first.chrs; i++)  {
    // on merge les bpoints 
    std::vector<double> pos(Z.first.bpoints[i].size() + Z.second.bpoints[i].size());
    std::merge(Z.first.bpoints[i].begin(), Z.first.bpoints[i].end(), 
               Z.second.bpoints[i].begin(), Z.second.bpoints[i].end(),
               pos.begin());

    double a = 0.;
    Z.first.set_cursor(i); Z.second.set_cursor(i);
    for(auto & b : pos) {
      Z.first.forward_cursor(b); Z.second.forward_cursor(b);
      if(HBD_at_point(Z))
        HBD += (b-a);
      a = b;
    }
  }
  return HBD;
}

}
