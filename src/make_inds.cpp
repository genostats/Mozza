//[[Rcpp::depends(gaston)]]
#include <Rcpp.h>
#include "mozza.h"
#include "gaston/matrix4.h"
using namespace Rcpp;

// cette fonction fait des individus indépendants
// avec des haplos mosaiques avec tuiles de longueurs 
// length_tiles parmi n_haps haplotypes
// Ils sont pusher dans un vecteur de mozza::zygote
void make_inds(std::vector<mozza::zygote> & ZYG, int n, int n_haps, double length_tiles = 20.) {
  for(int i = 0; i < n; i++) {
    ZYG.push_back( mozza::zygote(mozza::human_autosomes_b37, n_haps, length_tiles) );
  }
}

// idem avec l'autre constructeur de mosaiques
// qui prend un vecteur de proba des haplotypes de base
void make_inds(std::vector<mozza::zygote> & ZYG, int n, const std::vector<double> & proba_tiles, double length_tiles = 20.) {
  for(int i = 0; i < n; i++) {
    ZYG.push_back( mozza::zygote(mozza::human_autosomes_b37, proba_tiles, length_tiles) );
  }
}

// wrapers pour renvoyer des vecteurs de mozza::zygote
std::vector<mozza::zygote> make_inds(int n, int n_haps, double length_tiles = 20.) {
  std::vector<mozza::zygote> ZYG;
  make_inds(ZYG, n, n_haps, length_tiles);
  return ZYG;
}

std::vector<mozza::zygote> make_inds(int n, const std::vector<double> & proba_tiles, double length_tiles = 20.) {
  std::vector<mozza::zygote> ZYG;
  make_inds(ZYG, n, proba_tiles, length_tiles);
  return ZYG;
}


// Habillages avec un drop_to_bed_matrix pour finir
//[[Rcpp::export]]
List make_inds(int n, double length_tiles, XPtr<matrix4> Haplos, IntegerVector chr, NumericVector dist, bool kinship = false) {
  int n_haps = Haplos->ncol; // chaque haplotype = un "individu"
  std::vector<mozza::zygote> ZYG { make_inds(n, n_haps, length_tiles) }; 
  
  List L;
  L["bed"] = drop_to_bed_matrix(ZYG, Haplos, chr, dist);
  if(kinship) 
    L["kinship"] = kinship_matrix(ZYG);
  return L;
}


// proba_haplos est une matrice, chaque colonnes donne un jeu de proba sur les haplotypes
// le vecteur d'effectif N contient les effectifs à générer pour chacune des colonnees
// Cela permet de générer plus facilement des données avec des sous-populations ayant des 
// proportions différentes.
//[[Rcpp::export]]
List make_inds_probs(IntegerVector N, NumericMatrix proba_haplos, double length_tiles, 
                     XPtr<matrix4> Haplos, IntegerVector chr, NumericVector dist, bool kinship = false) {

  int n_haps = Haplos->ncol; // chaque haplotype = un "individu"
  if(n_haps != proba_haplos.nrow() || N.size() != proba_haplos.ncol())
    stop("Dimensions mismatch");

  std::vector<mozza::zygote> ZYG;
  
  std::vector<double> p( n_haps );
  for(int j = 0; j < proba_haplos.ncol(); j++) {
    for(int i = 0; i < n_haps; i++)
      p[i] = proba_haplos(i,j);
    make_inds(ZYG, N[j], p, length_tiles); 
  }

  List L;
  L["bed"] = drop_to_bed_matrix(ZYG, Haplos, chr, dist);
  if(kinship) 
    L["kinship"] = kinship_matrix(ZYG);
  return L;
}
