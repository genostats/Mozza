#include <Rcpp.h>
#include "mosaic.h"
#include "zygote.h"
#include "HBD_at_point.h"
#include "segments.h"
/* 
 * calcul des segments HBD
 */


namespace mozza {

// segments partagés HBD par les deux haplotypes du zygote
segments HBD_segments(zygote & Z) {
  segments HBD;
  for(int i = 0; i < Z.first.chrs; i++)  {
    // on merge les bpoints 
    std::vector<double> pos(Z.first.bpoints[i].size() + Z.second.bpoints[i].size());
    std::merge(Z.first.bpoints[i].begin(), Z.first.bpoints[i].end(), 
               Z.second.bpoints[i].begin(), Z.second.bpoints[i].end(),
               pos.begin());
    
    double a = 0.;
    Z.first.set_cursor(i); 
    Z.second.set_cursor(i);
    for(auto & b : pos) {
      Z.first.forward_cursor(b); 
      Z.second.forward_cursor(b);
      if(HBD_at_point(Z) && b > a) {
        HBD.chr.push_back(i);
        HBD.beg.push_back(a);
        HBD.end.push_back(b);
      }
      a = b;
    }
  }
  return HBD;
}

}
