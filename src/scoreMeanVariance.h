#include <Rcpp.h>
#include "gaston/matrix4.h"
#include "mapped_bed.h"

#ifndef _MOZZA_SCORE_EXPECTED_VALUE__
#define _MOZZA_SCORE_EXPECTED_VALUE__

namespace mozza {

// calcule la moyenne et la variance du score *pour un haplotype*
template<typename IV, typename DV>
std::pair<double,double> scoreMeanVariance(const mappedBed<IV, DV> & MB, const IV & submap, const DV & beta) {
  double m = 0, v = 0;
  int k = 0; // indice qui va être incrémenté au fur à mesure du parcours de la submap (cf beta[k])
  for(int i : submap) {
    // -- calcul du génotype moyen --
    unsigned int sumSNP = 0;
    int ncol = MB.haplotypes->ncol; // nb individus
    int true_ncol = MB.haplotypes->true_ncol;
    uint8_t * snp = MB.haplotypes->data[i];
    
    for(int ii = 0; ii < true_ncol-1; ii++) {
      unsigned char x = snp[ii];
      for(int ss = 0; ss < 4 && 4*ii+ss < ncol; ss++) {
        if(double(x&3) < 2) { // ce sont des haplotypes...
          sumSNP += (double) (x&3);
        }
        x >>= 2;
      }
    }
    double p = (double) sumSNP / ncol;
    // ---------
    double zz = beta[k] * p;
    m += zz;
    v += beta[k++] * zz * (1 -p); 
    
  }
  return std::make_pair(m, v); 
}
}
#endif