#include <Rcpp.h>
#include "mosaic.h"
#include "zygote.h"
namespace mozza {

// cette fonction n'a de sens que si les curseurs sont au même endroit dans les 4 haplotypes
// on ne vérifie **pas** ceci 
inline int IBD_at_point(zygote & Z1, zygote & Z2) {
  bool a = Z1.first.tile_at_cursor() == Z2.first.tile_at_cursor();
  bool b = Z1.first.tile_at_cursor() == Z2.second.tile_at_cursor();
  bool c = Z1.second.tile_at_cursor() == Z2.first.tile_at_cursor();
  bool d = Z1.second.tile_at_cursor() == Z2.second.tile_at_cursor();
  if(!a && !b && !c && !d)
    return 0;
  if((a && d) || (b && c))
    return 2;
  return 1;
}

// longueurs IBD 0, 1, et 2 entre Z1 et Z2
// note : on ne vérifie pas que la longueur des chromosomes / la fin de la dernière tuile coincident...
// NOTE la position des 4 curseurs est modifiée (on pourrait la sauver et remettre le curseur à sa place ?)
// NOTE pour rendre ça plus propre et lisible des fonctions qui bougent les deux curseurs d'un zygote en même
//      temps seraient bienvenues (voire des fonctions set_cursor(Z1, Z2, i), forward_cursor(Z1, Z2, b) ?)
std::tuple<double, double, double> IBD_length(zygote & Z1, zygote & Z2) {
  double IBD0 = 0, IBD1 = 0, IBD2 = 0.;
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
      switch(IBD_at_point(Z1, Z2)) {
      case 0:
        IBD0 += (b-a);
        break;
      case 1:
        IBD1 += (b-a);
        break;
      default:
        IBD2 += (b-a);
      }
      a = b;
    }
  }
  return std::make_tuple(IBD0, IBD1, IBD2);
}

}
