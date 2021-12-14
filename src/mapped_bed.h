#include <Rcpp.h>
#include "gaston/matrix4.h"

#ifndef _MAPPEDBED_
#define _MAPPEDBED_

namespace mozza {
template<typename IV, typename DB>
class mappedBed {
public:
  //   haplotypes = une bed matrix qui contient que des haplos (0, 1, ou 3 pour NA)
  //   chr  = vecteur de chromosomes (numérotés de 1 à ... : on leur soustrait 1 dans la fonction...)
  //   dist = un vecteur de positions (en cM)
  //   Ces deux vecteurs correspondent à la carte des SNPs présents dans haplotypes
  const XPtr<matrix4> haplotypes;
  const IV & chr; 
  const DB & dist;
  int nbSnps;

  mappedBed<IV, DB>(const XPtr<matrix4> Haps, const IV & chr_, const DB & dist_) : haplotypes(Haps), chr(chr_), dist(dist_) {
    nbSnps = chr.size();
    if(nbSnps != dist.size() || nbSnps != haplotypes->nrow)
      stop("Dimensions mismatch");
  }
};
}




#endif
