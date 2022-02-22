#include <Rcpp.h>
#include "gaston/matrix4.h"
#include "mapped_bed.h"

#ifndef _MOZZA_SCORE_ALL_HAPLOTYPES_
#define _MOZZA_SCORE_ALL_HAPLOTYPES_

namespace mozza {

// calcule le score de tous les haplotypes
// MB = mapped bed matrix (bed matrix + chr + dist)
// submap = vecteur d'indices des SNPs avec un effet
// beta = vecteur des effets
template<typename scalar_t, typename IV, typename DV>
std::vector<scalar_t> scoreAllHaplotypes(const mappedBed<IV, DV> & MB, const IV & submap, const DV & beta) {
  int ncol = MB.haplotypes->ncol; // nb individus = nb haplotypes
  int true_ncol = MB.haplotypes->true_ncol;

  std::vector<scalar_t> Scores(ncol);
  int k = 0; // indice de SNP qui va être incrémenté au fur à mesure du parcours de la submap (cf beta[k++])

  for(unsigned int i : submap) { // pour chaque SNP de la submap
    uint8_t * snp = MB.haplotypes->data[i];
    scalar_t bet = beta[k++];
    
    // on parcourt les individus / haplotypes 
    unsigned int in = 0; // indice d'individu
    for(unsigned int ii = 0; ii < true_ncol; ii++) {
      unsigned char x = snp[ii];
      for(int ss = 0; ss < 4 && 4*ii+ss < ncol; ss++) {
        if( (x&3) < 2) { // ce sont des haplotypes... on ne garde que les valeurs 0 ou 1
          Scores[ in++ ] += bet * (scalar_t) (x&3);
        }
        x >>= 2;
      }
    }
  }
  return Scores; 
}
}
#endif
