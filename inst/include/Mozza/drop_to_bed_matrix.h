//[[Rcpp::depends(gaston)]]
#include <Rcpp.h>
#include "gaston/matrix4.h"
#include "push_m4_row.h"
#include "mapped_bed.h"

using namespace Rcpp;

namespace mozza {

// x = un vecteur d'individus diploides (paires d'haplotypes)
// dans mapB il y a :
//   haplotypes = une bed matrix qui contient que des haplos (0, 1, ou 3 pour NA)
//   chr  = vecteur de chromosomes (numérotés de 1 à ... : on leur soustrait 1 dans la fonction...)
//   dist = un vecteur de positions (en cM)
//   Ces deux vecteurs correspondent à la carte des SNPs présents dans haplotypes
// On construit une bed matrix avec la même carte, qui contient les génotypes

template<typename IV, typename DB>
XPtr<matrix4> drop_to_bed_matrix(std::vector<zygote> & x, const mappedBed<IV, DB> mapB) {
  int nb_snps = mapB.nbSnps;
 
  int nb_inds = x.size(); // nb individus à créer
  XPtr<matrix4> A(new matrix4(nb_snps, nb_inds));
  
  // initialiser les positions
  int c = mapB.chr[0] - 1;
  for(auto & z : x) {
    z.first.set_cursor(c); 
    z.second.set_cursor(c);
  }
  // c'est parti
  for(int i = 0; i < nb_snps; i++) {
    // aller à la position du SNP (si besoin en changeant de chr)
    double pos = mapB.dist[i];
    if(c == mapB.chr[i] - 1) {
      for(auto & z : x) {
        z.first.forward_cursor(pos); 
        z.second.forward_cursor(pos);
      }
    } else {
      c = mapB.chr[i] - 1;
      for(auto & z : x) {
        z.first.set_cursor(c, pos); 
        z.second.set_cursor(c, pos);
      }   
    }
    // alleles pour le SNP dans haplotypes->data[i]
    // pour pas s'embêter on va juste les copier dans un std::vector
    std::vector<char> alleles;
    push_m4_row(mapB.haplotypes->data[i], mapB.haplotypes->ncol, alleles);
    // l'objet S permet de push back des valeurs dans A->data[i]
    SNP_push_back S(A->data[i], A->ncol);
    // la fonction push back 'dans S' les génotypes à la position courante du curseur
    push_genotypes_at_cursor(x, alleles, S);
  }
  return A;
}


// Idem avec seulement des haplotypes
template<typename IV, typename DB>
XPtr<matrix4> drop_to_bed_matrix(std::vector<mosaic> & x, const mappedBed<IV, DB> mapB) {
  int nb_snps = mapB.nbSnps;
  int nb_inds = x.size(); // nb "individus" (haplotypes) à créer
  XPtr<matrix4> A(new matrix4(nb_snps, nb_inds));
  
  // initialiser les positions
  int c = mapB.chr[0] - 1;
  for(auto & z : x) {
    z.set_cursor(c); 
  }
  // c'est parti
  for(int i = 0; i < nb_snps; i++) {
    // aller à la position du SNP (si besoin en changeant de chr)
    double pos = mapB.dist[i];
    if(c == mapB.chr[i] - 1) {
      for(auto & z : x) {
        z.forward_cursor(pos); 
      }
    } else {
      c = mapB.chr[i] - 1;
      for(auto & z : x) {
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
    push_genotypes_at_cursor(x, alleles, S);
  }
  return A;
}

}
