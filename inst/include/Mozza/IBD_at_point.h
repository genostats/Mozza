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

}

