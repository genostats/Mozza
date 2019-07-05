//[[Rcpp::depends(gaston)]]
#include <Rcpp.h>
#include "mozza.h"
using namespace Rcpp;

//[[Rcpp::export]]
List nuclear_families(int nb_fams, int nb_offsprings, double tile_length, XPtr<matrix4> Haplos, IntegerVector chr, NumericVector dist, bool kinship = false) {
  std::vector<mozza::zygote> ZYG; 
  int n_haps = Haplos->ncol; // chaque haplotype = un "individu"
  
  for(int i = 0; i < nb_fams; i++) {
    mozza::zygote M(mozza::human_autosomes_b37, n_haps, tile_length);
    mozza::zygote F(mozza::human_autosomes_b37, n_haps, tile_length);
    ZYG.push_back(M);
    ZYG.push_back(F);
    for(int j = 0; j < nb_offsprings; j++)
      ZYG.push_back(M+F);
  }
  List L;
  L["bed"] = drop_to_bed_matrix(ZYG, Haplos, chr, dist);
  if(kinship) 
    L["kinship"] = kinship_matrix(ZYG);
  return L;
}


// Version avec haplotypes probabilisés
// nb_fams familles à nb_offsprings
void nuclear_families(std::vector<mozza::zygote> & ZYG, int nb_fams, int nb_offsprings, const std::vector<double> & proba_tiles,
                      double tile_length) {
  for(int i = 0; i < nb_fams; i++) {
    mozza::zygote M(mozza::human_autosomes_b37, proba_tiles, tile_length);
    mozza::zygote F(mozza::human_autosomes_b37, proba_tiles, tile_length);
    ZYG.push_back(M);
    ZYG.push_back(F);
    for(int j = 0; j < nb_offsprings; j++)
      ZYG.push_back(M+F);
  }
}

// ce code ressemble à celui de make_inds_probs (voir commentaires à cette fction)
// Nfams = autant d'élts que de colonnes à proba_haplos
// chaque famille a nb_offsprings
//[[Rcpp::export]]
List nuclear_families_probs(IntegerVector Nfams, int nb_offsprings, NumericMatrix proba_haplos, double length_tiles,
                     XPtr<matrix4> Haplos, IntegerVector chr, NumericVector dist, bool kinship = false) {

  int n_haps = Haplos->ncol; // chaque haplotype = un "individu"
  if(n_haps != proba_haplos.nrow() || Nfams.size() != proba_haplos.ncol())
    stop("Dimensions mismatch");

  std::vector<mozza::zygote> ZYG;

  std::vector<double> p( n_haps );
  for(int j = 0; j < proba_haplos.ncol(); j++) {
    for(int i = 0; i < n_haps; i++)
      p[i] = proba_haplos(i,j);
    nuclear_families(ZYG, Nfams[j], nb_offsprings, p, length_tiles);
  }

  List L;
  L["bed"] = drop_to_bed_matrix(ZYG, Haplos, chr, dist);
  if(kinship)
    L["kinship"] = kinship_matrix(ZYG);
  return L;
}

