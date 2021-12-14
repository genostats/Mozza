//[[Rcpp::depends(gaston)]]
#include <Rcpp.h>
#include "gaston/matrix4.h"
#include "mozza.h"
#include "getHaploScore.h"

using namespace Rcpp;

#ifndef __GETZYGOTESCORE__
#define __GETZYGOTESCORE__

namespace mozza {

template<typename INV1, typename DBV1, typename INV2, typename DBV2>
double getZygoteScore(zygote & z, XPtr<matrix4> haplotypes, const INV1 & chr, const DBV1 & dist,
                     const INV2 & submap, const DBV2 & beta) {
  return getHaploScore(z.first, haplotypes, chr, dist, submap, beta) + getHaploScore(z.second, haplotypes, chr, dist, submap, beta);
}

}


#endif

