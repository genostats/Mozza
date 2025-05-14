#include <Rcpp.h>
#include "gaston/matrix4.h"

#ifndef _MAPPEDBED_
#define _MAPPEDBED_

// pour une raison que j'ai oubliée, je ne peux pas utiliser DV comme typename (?? est ce bien vrai ?)
namespace mozza {
template<typename IV, typename DB, typename CV = std::vector<std::string>>
class mappedBed {
public:
  //   haplotypes = une bed matrix qui contient que des haplos (0, 1, ou 3 pour NA)
  //   chr  = vecteur de chromosomes (numérotés de 1 à ... : on leur soustrait 1 dans la fonction...)
  //   dist = un vecteur de positions (en cM)
  //   Ces deux vecteurs correspondent à la carte des SNPs présents dans les haplotypes
  const XPtr<matrix4> haplotypes;
  const IV & chr; 
  const DB & dist;
  
  // membres facultatifs [pour le drop to vcf...]
  const CV & id = {};
  const IV & pos = {};
  const CV & A1 = {};
  const CV & A2 = {};

  unsigned int nbSnps;
  unsigned int nbHaps;

  mappedBed<IV, DB, CV>(const XPtr<matrix4> Haps, const IV & chr_, const DB & dist_) : haplotypes(Haps), chr(chr_), dist(dist_) {
    nbSnps = chr.size();
    nbHaps = Haps->ncol;  
    if(nbSnps != dist.size() || nbSnps != haplotypes->nrow)
      stop("Dimensions mismatch");
  }

  mappedBed<IV, DB, CV>(const XPtr<matrix4> Haps, const IV & chr_, const DB & dist_, const CV & id_, const IV & pos_,
                        const CV & A1_, const CV & A2_) : 
                        haplotypes(Haps), chr(chr_), dist(dist_), id(id_), pos(pos_), A1(A1_), A2(A2_)
  {
    nbSnps = chr.size();
    nbHaps = Haps->ncol;  
    if(nbSnps != dist.size() || nbSnps != haplotypes->nrow)
      stop("Dimensions mismatch");
  }


};
}




#endif
