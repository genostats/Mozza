#ifndef MOZZA_
#define MOZZA_

#include "gaston/matrix4.h"
#include "Mozza/mosaic.h"
#include "Mozza/human_autosomes.h"
#include "Mozza/push_genotypes_at_cursor.h"
#include "Mozza/push_haplotypes_at_cursor.h"
#include "Mozza/zygote.h"
#include "Mozza/mapped_bed.h"
#include "Mozza/zygote_to_bed_matrix.h"
#include "Mozza/haplo_to_bed_matrix.h"
#include "Mozza/kinship_matrix.h"
#include "Mozza/Autosomes.h"

namespace mozza {

  double IBD_sharing(mosaic & M1, mosaic & M2);

  inline int IBD_at_point(zygote & Z1, zygote & Z2);
  std::tuple<double, double, double> IBD_length(zygote & Z1, zygote & Z2);

  inline bool HBD_at_point(zygote & Z);
  double HBD_length(zygote & Z);

  inline bool relatedness_at_point(zygote & Z);
  double relatednessLength(zygote &, zygote &);

  NumericMatrix ibd_matrix(std::vector<mosaic> & HAP);

//  NumericMatrix kinship_matrix(std::vector<zygote> & ZYG);
  NumericMatrix fraternity_matrix(std::vector<zygote> & ZYG);
}

#endif

