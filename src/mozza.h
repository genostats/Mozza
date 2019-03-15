#ifndef MOZZA_
#define MOZZA_

#include "mosaic.h"
#include "human_autosomes.h"
#include "push_genotypes_at_cursor.h"
#include "gaston/matrix4.h"
#include "zygote.h"

namespace mozza {

  double sharing(mosaic & M1, mosaic & M2);
  XPtr<matrix4> drop_to_bed_matrix(std::vector<std::pair<mosaic,mosaic>> & x, XPtr<matrix4> haplotypes, IntegerVector chr, NumericVector dist);

  inline int IBD_at_point(zygote & Z1, zygote & Z2);
  std::tuple<double, double, double> IBD_length(zygote & Z1, zygote & Z2);
    
};

#endif