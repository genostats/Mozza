#ifndef MOZZA_
#define MOZZA_

#include "mosaic.h"
#include "human_autosomes.h"
#include "push_genotypes_at_cursor.h"
#include "gaston/matrix4.h"
#include "zygote.h"
#include "mapped_bed.h"
#include "drop_to_bed_matrix.h"

namespace mozza {

  double sharing(mosaic & M1, mosaic & M2);
  
  inline int IBD_at_point(zygote & Z1, zygote & Z2);
  std::tuple<double, double, double> IBD_length(zygote & Z1, zygote & Z2);
    
  inline bool HBD_at_point(zygote & Z);
  double HBD_length(zygote & Z);

  inline bool relatedness_at_point(zygote & Z);
  double relatednessLength(zygote &, zygote &);
     
  NumericMatrix ibd_matrix(std::vector<mosaic> & HAP);

  NumericMatrix kinship_matrix(std::vector<zygote> & ZYG);
  NumericMatrix fraternity_matrix(std::vector<zygote> & ZYG);
}

#endif
