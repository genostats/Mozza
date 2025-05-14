#include <Rcpp.h>
#include "mosaic.h"
#include "zygote.h"
#include "getZygote.h"

using namespace Rcpp;

#ifndef PUSH_GENOTYPES_AT_CURSOR_
#define PUSH_GENOTYPES_AT_CURSOR_

namespace mozza {
// on fait un template car comme ça on pourra passer
// n'importe quoi qui a une méthode push_back() pour ranger les génotypes
// donc par exemple un std::vector 
// mais on peut fabriquer une classe avec une méthode push_back pour 
// remplir directement un char * dans une bed matrix
// idem pour la lecture des génotypes, n'importe quoi qui a une méthode []
// fera l'affaire
// ça suppose plus ou moins implicitement que les curseur sont à la même
// position partout mais peu importe
template<typename ZV, typename T1, typename T2>
inline void push_genotypes_at_cursor(ZV & x, T1 & allele_at_haplo, T2 & genotypes) {
  for(auto & pp : x) {
    zygote & p(getZygote(pp));
    int h1 = p.first.tile_at_cursor();
    int h2 = p.second.tile_at_cursor();
    genotypes.push_back( allele_at_haplo[h1], allele_at_haplo[h2] );
  }
}

}

#endif
