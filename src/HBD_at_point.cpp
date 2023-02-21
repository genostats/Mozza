#include <Rcpp.h>
#include "mosaic.h"
#include "zygote.h"

namespace mozza {

// cette fonction n'a de sens que si les curseurs sont au même endroit dans les 2 haplotypes
// on ne vérifie **pas** ceci 
bool HBD_at_point(zygote & Z) {
  return(Z.first.tile_at_cursor() == Z.second.tile_at_cursor());
}

}
