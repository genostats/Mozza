//[[Rcpp::depends(gaston)]]
#include "Mozza/mosaic.h"
#include "Mozza/zygote.h"
#include "gaston/matrix4.h"
#include "push_m4_row.h"
#include "mapped_bed.h"
#include "getHaplo.h"

#ifndef _mozza_haplo_to_bed_matrix_
#define _mozza_haplo_to_bed_matrix_

namespace mozza {

template<typename HV, typename IV, typename DB>
XPtr<matrix4> haplo_to_bed_matrix(HV & HAP, const mappedBed<IV, DB> mapB) {
  int nb_snps = mapB.nbSnps;
  int nb_inds = HAP.size(); // nb "individus" (haplotypes) à créer
  XPtr<matrix4> A(new matrix4(nb_snps, nb_inds));
  
  // initialiser les positions
  int c = mapB.chr[0] - 1;
  for(int k = 0; k < nb_inds; k++) {
    mosaic & z(getHaplo(HAP[k]));
    z.set_cursor(c); 
  }
  // c'est parti
  for(int i = 0; i < nb_snps; i++) {
    // aller à la position du SNP (si besoin en changeant de chr)
    double pos = mapB.dist[i];
    if(c == mapB.chr[i] - 1) { // pas de chgt de chr : forward cursor
      for(int k = 0; k < nb_inds; k++) {
        mosaic & z(getHaplo(HAP[k]));
        z.forward_cursor(pos); 
      }
    } else {
      c = mapB.chr[i] - 1;
      for(int k = 0; k < nb_inds; k++) { // nouveau chr : set cursor
        mosaic & z(getHaplo(HAP[k]));
        z.set_cursor(c, pos); 
      }   
    }
    // alleles pour le SNP dans haplotypes->data[i]
    // pour pas s'embêter on va juste les copier dans un std::vector
    std::vector<char> alleles;
    push_m4_row(mapB.haplotypes->data[i], mapB.haplotypes->ncol, alleles);
    // l'objet S permet de push back des valeurs dans A->data[i]
    SNP_push_back S(A->data[i], A->ncol);
    // la fonction push back dans S les alleles à la position courante du curseur
    push_haplotypes_at_cursor(HAP, alleles, S);
  }
  return A;
}

}
#endif
