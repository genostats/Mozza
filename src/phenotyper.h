#include <Rcpp.h>
#include "gaston/matrix4.h"
#include "mapped_bed.h"
#include "scoreAllHaplotypes.h"

#ifndef _MOZZA_PHENOTYPER__
#define _MOZZA_PHENOTYPER_

namespace mozza {

template<typename IV, typename DV>
class phenotyper {
  const mappedBed<IV, DV> & MB; // (bed matrix + chr + dist)
  const IV & submap;            //  vecteur d'indices des SNPs avec un effet
  const DV & beta;              //  vecteur des effets
public:
  double mu;                    //  moyenne du score zygotique
  double s, sdE;                //  s = pour mettre à la bonne échelle les "raw scores", sdE = variance de l'environnement
  // constructeur 1, utilise scoreAllHaplotypes pour avoir une estimation de l'espérance et de la variance des scores
  // pas évident que ça soit pas un peu surestimé puisqu'en pratique on va brasser les haplotypes ! 
  // L'autre solution est d'utiliser scoreMeanVariance qui fait le même calcul en supposant l'équilibre de liaison...

  phenotyper<IV, DV>(const mappedBed<IV, DV> & MB_, const IV & submap_, const DV & beta_) 
                     : MB(MB_), submap(submap_), beta(beta_) {
    if(beta.size() != submap.size()) {
      stop("Submap size and beta coeff size mismatch");
    }
    // uncalibrated..
    mu = 0; s = 1; sdE = 1;
  }

  void setCalibration(double mu_, double s_, double sdE_) {
    mu = mu_;
    s = s_;
    sdE = sdE_;
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
      score += (double) MB.haplotypes->get(i, x.tile_at_cursor()) * beta[k++];
    }
    return score;
  }

  // cette fonction renvoie un score "brut"
  double getZygoteScore(mozza::zygote & z) {
    return getHaploScore(z.first) + getHaploScore(z.second);
  }

  // et celle ci met à l'échelle pour que la liabilité soit centrée réduite
  std::tuple<double, double, double> getLiability(mozza::zygote & z) {
    double G0 = getZygoteScore(z);
    double G1 = (G0 - mu)*s;
    double E = R::norm_rand()*sdE;
    return std::make_tuple(G0, G1, E);
  }

};

}
#endif
