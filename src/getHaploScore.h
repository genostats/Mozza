//[[Rcpp::depends(gaston)]]
#include <Rcpp.h>
#include "gaston/matrix4.h"
#include "mozza.h"

using namespace Rcpp;

#ifndef __GETHAPLOSCORE__
#define __GETHAPLOSCORE__
namespace mozza {

// x = un vecteur d'individus diploides (paires d'haplotypes)
// haplotypes = une bed matrix qui contient que des haplos (0, 1, ou 3 pour NA)
// chr  = vecteur de chromosomes (numérotés de 1 à ... : on leur soustrait 1 dans la fonction...)
// dist = un vecteur de positions (en cM)
// Ces deux vecteurs correspondent à la carte des SNPs présents dans haplotypes
// submap = un vecteur d'incides (entre 0 et nb_snps - 1)
// beta = vecteur d'effets alléliques
template<typename INV1, typename DBV1, typename INV2, typename DBV2>
double getHaploScore(mosaic & x, XPtr<matrix4> haplotypes, const INV1 & chr, const DBV1 & dist, 
                     const INV2 & submap, const DBV2 & beta) {

  int nb_snps_total = chr.size();
  if(nb_snps_total != dist.size() || nb_snps_total != haplotypes->nrow)
    stop("Map and bed matrix dimensions mismatch");
  
  int nb_snps = submap.size();
  if(nb_snps != beta.size())
    stop("Submap size and beta coeff size mismatch");

  // initialiser les positions
  int c = chr[0] - 1;
  x.set_cursor(c);
  
  double score = 0;
  int k = 0; // il faut parcourir à la fois les élts de submap et ceux de beta
             // on a opté pour cette solution (cf derniere ligne de la boucle)
  // c'est parti
  for(int i : submap) {
    // aller à la position du SNP (si besoin en changeant de chr)
    double pos = dist[i];
    if(c == chr[i] - 1) {
      x.forward_cursor(pos);
    } else { // changer de chromosome
      c = chr[i] - 1;
      x.set_cursor(c, pos);
    }
    // les alleles pour le SNP sont dans haplotypes->data[i]
    // contrairement à la fonction 'drop_to_bed_matrix' on n'a besoin
    // que d'un seul d'entre eux !
    score += (double) haplotypes->get(i, x.tile_at_cursor()) * beta[k++];
  }
  return score;
}

}
#endif

