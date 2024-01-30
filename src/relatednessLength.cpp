#include <Rcpp.h>
#include "mosaic.h"
#include "zygote.h"
/* 
 *
 * calcul de la relatedness moyenne
 * écrit à partir des fonctions IBD_at_point et IBD_length (de ibd.cpp)
 *
 */


namespace mozza {

// cette fonction n'a de sens que si les curseurs sont au même endroit dans les 4 haplotypes
// on ne vérifie **pas** ceci 
inline double relatedness_at_point(zygote & Z1, zygote & Z2) {
  int a = Z1.first.tile_at_cursor();
  int b = Z1.second.tile_at_cursor();
  int c = Z2.first.tile_at_cursor();
  int d = Z2.second.tile_at_cursor();
  int x = (int) (a == c) + (int) (a == d) + (int) (b == c) + (int) (b == d);
  return 0.25 * (double) x;
}

// somme des relatedness x longueur des segments...
// note : on ne vérifie pas que la longueur des chromosomes / la fin de la dernière tuile coincident...
//      cf ibd.cpp pour plus de réflexions sur la façon de rendre ça propre
double relatednessLength(zygote & Z1, zygote & Z2) {
  double RL = 0;
  for(int i = 0; i < Z1.first.chrs; i++)  {
    // on merge les bpoints en trois étapes
    // possibilité de faire mieux ?
    std::vector<double> pos1(Z1.first.bpoints[i].size() + Z1.second.bpoints[i].size());
    std::vector<double> pos2(Z2.first.bpoints[i].size() + Z2.second.bpoints[i].size());
    std::vector<double> pos(pos1.size() + pos2.size());
    std::merge(Z1.first.bpoints[i].begin(), Z1.first.bpoints[i].end(), 
               Z1.second.bpoints[i].begin(), Z1.second.bpoints[i].end(),
               pos1.begin());

    std::merge(Z2.first.bpoints[i].begin(), Z2.first.bpoints[i].end(), 
               Z2.second.bpoints[i].begin(), Z2.second.bpoints[i].end(),
               pos2.begin());
    
    std::merge(pos1.begin(), pos1.end(), pos2.begin(), pos2.end(), pos.begin());
    
    double a = 0.;
    Z1.first.set_cursor(i); Z1.second.set_cursor(i);
    Z2.first.set_cursor(i); Z2.second.set_cursor(i);
    for(auto & b : pos) {
      Z1.first.forward_cursor(b); Z1.second.forward_cursor(b);
      Z2.first.forward_cursor(b); Z2.second.forward_cursor(b);
      
      RL += (b-a)*relatedness_at_point(Z1, Z2);
      a = b;
    }
  }
  return RL;
}

}
