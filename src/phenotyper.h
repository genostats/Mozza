#include <Rcpp.h>
#include "gaston/matrix4.h"
#include "mapped_bed.h"

#ifndef _MOZZA_PHENOTYPER__
#define _MOZZA_PHENOTYPER_

namespace mozza {

template<typename IV, typename DV>
class phenotyper {
  const mappedBed<IV, DV> & MB;
  const IV & submap;
  const DV & beta;
public:
  phenotyper<IV, DV>(const mappedBed<IV, DV> & MB_, const IV & submap_, const DV & beta_) : 
  MB(MB_), submap(submap_), beta(beta_) {
    if(beta.size() != submap.size()) {
      stop("Submap size and beta coeff size mismatch");
    }
  }

  double getHaploScore(mozza::mosaic & x) {
    double score = 0;

    // initialiser les positions
    int c = MB.chr[0] - 1;
    x.set_cursor(c);

    int k = 0; // il faut parcourir à la fois les élts de submap et ceux de beta
               // on a opté pour la solution qui consiste à parcourir submap en 
               // incrémentant l'indice k à chaque tour (cf derniere ligne de la boucle)
    // c'est parti
    for(int i : submap) {
      // aller à la position du SNP (si besoin en changeant de chr)
      double pos = MB.dist[i];
      if(c == MB.chr[i] - 1) {
        x.forward_cursor(pos);
      } else { // changer de chromosome
        c = MB.chr[i] - 1;
        x.set_cursor(c, pos);
      }
      // les alleles pour le SNP sont dans haplotypes->data[i]
      // contrairement à la fonction 'drop_to_bed_matrix' on n'a besoin
      // que d'un seul d'entre eux !
// Rcpp::Rcout << (int) MB.haplotypes->get(i, x.tile_at_cursor()) << "x" << beta[k] << " + ";
      score += (double) MB.haplotypes->get(i, x.tile_at_cursor()) * beta[k++];
    }
    return score;
  }

  double getZygoteScore(mozza::zygote & z) {
    return getHaploScore(z.first) + getHaploScore(z.second);
  }
};

}
#endif
