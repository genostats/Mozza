#include <Rcpp.h>
#include "gaston/matrix4.h"
#include "mapped_bed.h"

#ifndef _MOZZA_SCORE_EXPECTED_VALUE__
#define _MOZZA_SCORE_EXPECTED_VALUE__

namespace mozza {

// calcule la moyenne et la variance du score *pour un haplotype*

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!! la variance est calculée en supposant l'équilibre de liaison !!!!!!
// !!!!!! ça n'est pas une (très) bonne approximation                  !!!!!!
// !!!!!! on garde ce morceau de code pour le moment mais on va plutôt !!!!!!
// !!!!!! utiliser scoreAllHaplotypes pour avoir un truc exact         !!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

// MB = mapped bed matrix (bed matrix + chr + dist)
// submap = vecteur d'indices des SNPs avec un effet
// beta = vecteur des effets
template<typename IV, typename DV>
std::pair<double,double> scoreMeanVariance(const mappedBed<IV, DV> & MB, const IV & submap, const DV & beta) {
  int ncol = MB.haplotypes->ncol; // nb individus
  int true_ncol = MB.haplotypes->true_ncol;

  double m = 0, v = 0;
  int k = 0; // indice qui va être incrémenté au fur à mesure du parcours de la submap (cf beta[k])

  for(int i : submap) { // pour chaque SNP de la submap
    // -- calcul du génotype moyen --
    uint8_t * snp = MB.haplotypes->data[i];
    
    // on calcule p = freq de l'allèle alt 
    unsigned int sumSNP = 0; 
    for(int ii = 0; ii < true_ncol; ii++) {
      unsigned char x = snp[ii];
      for(int ss = 0; ss < 4 && 4*ii+ss < ncol; ss++) {
        if( (x&3) < 2) { // ce sont des haplotypes... on ne garde que les valeurs 0 ou 1
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
