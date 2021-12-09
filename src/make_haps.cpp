//[[Rcpp::depends(gaston)]]
#include <Rcpp.h>
#include "mozza.h"
#include "gaston/matrix4.h"
using namespace Rcpp;

// cette fonction fait des individus indépendants
// avec des haplos mosaiques avec tuiles de longueurs 
// length_tiles parmi n_haps haplotypes
// Ils sont pushed back dans un vecteur de mozza::mosaic
void make_haps(std::vector<mozza::mosaic> & HAP, int n, int n_haps, double length_tiles = 20.) {
  for(int i = 0; i < n; i++) {
    HAP.push_back( mozza::mosaic(mozza::human_autosomes_b37, n_haps, length_tiles) );
  }
}

// idem avec l'autre constructeur de mosaiques
// qui prend un vecteur de proba des haplotypes de base
void make_haps(std::vector<mozza::mosaic> & HAP, int n, const std::vector<double> & proba_tiles, double length_tiles = 20.) {
  for(int i = 0; i < n; i++) {
    HAP.push_back( mozza::mosaic(mozza::human_autosomes_b37, proba_tiles, length_tiles) );
  }
}

// wrapers pour renvoyer des vecteurs de mozza::mosaic
std::vector<mozza::mosaic> make_haps(int n, int n_haps, double length_tiles = 20.) {
  std::vector<mozza::mosaic> HAP;
  make_haps(HAP, n, n_haps, length_tiles);
  return HAP;
}

std::vector<mozza::mosaic> make_haps(int n, const std::vector<double> & proba_tiles, double length_tiles = 20.) {
  std::vector<mozza::mosaic> HAP;
  make_haps(HAP, n, proba_tiles, length_tiles);
  return HAP;
}


// Habillages avec un drop_to_bed_matrix pour finir
//[[Rcpp::export]]
List make_haps(int n, double length_tiles, XPtr<matrix4> Haplos, IntegerVector chr, NumericVector dist, 
               bool ibd = false) {
  int n_haps = Haplos->ncol; // chaque haplotype = un "individu" de la bed matrix
  std::vector<mozza::mosaic> HAP { make_haps(n, n_haps, length_tiles) }; 
  
  List L;
  L["bed"] = drop_to_bed_matrix(HAP, Haplos, chr, dist);
  if(ibd) 
    L["ibd"] = ibd_matrix(HAP);
  return L;
}


// proba_haplos est une matrice, chaque colonnes donne un jeu de proba sur les haplotypes
// le vecteur d'effectif N contient les effectifs à générer pour chacune des colonnees
// Cela permet de générer plus facilement des données avec des sous-populations ayant des 
// proportions différentes.
//[[Rcpp::export]]
List make_haps_probs(IntegerVector N, NumericMatrix proba_haplos, double length_tiles, 
                     XPtr<matrix4> Haplos, IntegerVector chr, NumericVector dist, 
                     bool ibd = false) {

  int n_haps = Haplos->ncol; // chaque haplotype = un "individu" de la bed matrix
  if(n_haps != proba_haplos.nrow() || N.size() != proba_haplos.ncol()) {
    Rcerr << "nb haplotypes = " << n_haps << "\n";
    Rcerr << "probability matrix has " << proba_haplos.nrow() << " rows\n";
    Rcerr << "probability matrix has " << proba_haplos.ncol() << " cols\n";
    Rcerr << "number of demes " << N.size() << "\n";
    stop("Dimensions mismatch");
  }
  std::vector<mozza::mosaic> HAP;
  
  std::vector<double> p( n_haps );
  for(int j = 0; j < proba_haplos.ncol(); j++) {
    for(int i = 0; i < n_haps; i++)
      p[i] = proba_haplos(i,j);
    make_haps(HAP, N[j], p, length_tiles); 
  }

  List L;
  L["bed"] = drop_to_bed_matrix(HAP, Haplos, chr, dist);
  if(ibd)
    L["ibd"] = ibd_matrix(HAP);
  return L;
}
