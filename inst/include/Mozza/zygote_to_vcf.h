//[[Rcpp::depends(gaston)]]
#include "Mozza/mosaic.h"
#include "Mozza/zygote.h"
#include "push_m4_row.h"
#include "mapped_bed.h"
#include "getZygote.h"
#include <fstream>

#ifndef _mozza_zygote_to_vcf_
#define _mozza_zygote_to_vcf_

namespace mozza {

// x = un vecteur d'individus diploides (paires d'haplotypes)
// dans mapB il y a :
//   haplotypes = une bed matrix qui contient que des haplos (0, 1, ou 3 pour NA)
//   chr  = vecteur de chromosomes (numérotés de 1 à ... : on leur soustrait 1 dans la fonction...)
//   dist = un vecteur de positions (en cM)
//   Ces deux vecteurs correspondent à la carte des SNPs présents dans haplotypes
// On écrit les génotypes phasés dans un vcf

// !!!!! Cette fonction est quasi identique à zygote_to_bed_matrix : peut-on refactoriser ce code ?

// ZV = Zygote Vector (vecteur de zygotes *ou* de pointeurs grâce à getZygote)
// IV = Integer Vector
// DB = Double Vector
// CV = Character Vector
template<typename ZV, typename IV, typename DB, typename CV>
void zygote_to_vcf(ZV & ZYG, const mappedBed<IV, DB, CV> & mapB, std::ofstream & output) {
  int nb_snps = mapB.nbSnps;
 
  int nb_inds = ZYG.size(); // nb individus à créer

  // initialiser les positions
  int c = mapB.chr[0] - 1;
  for(int k = 0; k < nb_inds; k++) {
    zygote & z(getZygote(ZYG[k]));
    z.first.set_cursor(c);
    z.second.set_cursor(c);
  }
  // c'est parti
  for(int i = 0; i < nb_snps; i++) {
    // aller à la position du SNP (si besoin en changeant de chr)
    double pos = mapB.dist[i];
    if(c == mapB.chr[i] - 1) { // on est sur le bon chr : forward cursor
      for(int k = 0; k < nb_inds; k++) {
        zygote & z(getZygote(ZYG[k]));
        z.first.forward_cursor(pos);
        z.second.forward_cursor(pos);
      }
    } else { // on change de chr (set cursor)
      c = mapB.chr[i] - 1;
      for(int k = 0; k < nb_inds; k++) {
        zygote & z(getZygote(ZYG[k]));
        z.first.set_cursor(c, pos);
        z.second.set_cursor(c, pos);
      }
    }
    // alleles pour le SNP dans haplotypes->data[i]
    // pour pas s'embêter on va juste les copier dans un std::vector
    std::vector<char> alleles;
    push_m4_row(mapB.haplotypes->data[i], mapB.haplotypes->ncol, alleles);
    // l'objet S permet de push back des valeurs dans le fichier, au format VCF
    SNP_push_back_to_vcf S(output, mapB.chr[i], mapB.pos[i], as<std::string>(mapB.id[i]), as<std::string>(mapB.A1[i]), as<std::string>(mapB.A2[i]));
    // la fonction push back 'dans S' les génotypes à la position courante du curseur
    push_genotypes_at_cursor(ZYG, alleles, S);
  }
}

}

#endif
