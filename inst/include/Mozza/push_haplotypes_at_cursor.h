#include <Rcpp.h>
#include "mosaic.h"
#include "zygote.h"
#include "getHaplo.h"

using namespace Rcpp;

#ifndef PUSH_HAPLOTYPES_AT_CURSOR_
#define PUSH_HAPLOTYPES_AT_CURSOR_
// cf commentaires dans push_genotypes_at_cursor.h

namespace mozza {
template<typename HV, typename T1, typename T2>
inline void push_haplotypes_at_cursor(HV & x, T1 & allele_at_haplo, T2 & genotypes) {
  for(auto & pp : x) {
    mosaic & p(getHaplo(pp));
    int h = p.tile_at_cursor();
    genotypes.push_back( allele_at_haplo[h] );
  }
}
}

#endif
