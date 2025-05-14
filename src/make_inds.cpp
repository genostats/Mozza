//[[Rcpp::depends(gaston)]]
#include <Rcpp.h>
#include "Mozza.h"
#include "gaston/matrix4.h"
using namespace Rcpp;

// cette fonction fait des individus indépendants
// avec des haplos mosaiques avec tuiles de longueurs 
// length_tiles parmi n_haps haplotypes
// Ils sont pushed back dans un vecteur de mozza::zygote
// ZYG : vector in which to push pack the zygote
// n : number of inds to generate
// n_haps : tiles will be numbered from 0 to n_haps-1
// length_tiles : in cM
void make_inds(std::vector<mozza::zygote> & ZYG, int n, int n_haps, double length_tiles = 20.) {
  for(int i = 0; i < n; i++) {
    ZYG.push_back( mozza::zygote(mozza::Autosomes(), n_haps, length_tiles) );
  }
}

// idem avec un vecteur de pointeurs
void make_inds(std::vector<mozza::zygote *> & ZYG, int n, int n_haps, double length_tiles = 20.) {
  for(int i = 0; i < n; i++) {
    mozza::zygote * pz = new mozza::zygote(mozza::Autosomes(), n_haps, length_tiles);
    ZYG.push_back( pz );
  }
}

// R export
//[[Rcpp::export]]
List make_inds(int n, unsigned int n_haps, double length_tiles) {
  std::vector<mozza::zygote *> x; 
  make_inds(x, n, n_haps, length_tiles);
  
  List ZYG(n);
  for(int i = 0; i < n; i++) ZYG[i] = XPtr<mozza::zygote>(x[i], true);
  ZYG.attr("class") = CharacterVector::create("zygote");
  return ZYG;
}

// --------------------------------------------------------------------------------------------

// idem avec l'autre constructeur de mosaiques
// qui prend un vecteur de proba des haplotypes de base
// ZYG : vector in which to push pack the zygote
// n : number of inds to generate
// proba_tiles : vector of probabilities (its length n_hap give the bound of the tile numbers)
// length_tiles : in cM
void make_inds(std::vector<mozza::zygote> & ZYG, int n, const std::vector<double> & proba_tiles, double length_tiles = 20.) {
  for(int i = 0; i < n; i++) {
    ZYG.push_back( mozza::zygote(mozza::Autosomes(), proba_tiles, length_tiles) );
  }
}

// idem avec un vecteur de pointeurs
void make_inds(std::vector<mozza::zygote *> & ZYG, int n, const std::vector<double> & proba_tiles, double length_tiles = 20.) {
  for(int i = 0; i < n; i++) {
    mozza::zygote * pz = new mozza::zygote(mozza::Autosomes(), proba_tiles, length_tiles);
    ZYG.push_back( pz );
  }
}

// proba_haplos est une matrice, chaque colonnes donne un jeu de proba sur les haplotypes
// le vecteur d'effectif N contient les effectifs à générer pour chacune des colonnes
// Cela permet de générer plus facilement des données avec des sous-populations ayant des 
// proportions différentes.
//[[Rcpp::export]]
List make_inds_probs(IntegerVector N, NumericMatrix proba_haplos, double length_tiles) {

  int n_haps = proba_haplos.nrow(); // chaque haplotype = un "individu"
  if(N.size() != proba_haplos.ncol()) {
    Rcerr << "probability matrix has " << proba_haplos.ncol() << " cols\n";
    Rcerr << "number of demes " << N.size() << "\n";
    stop("Dimensions mismatch");
  }
  std::vector<mozza::zygote *> x;
  
  std::vector<double> p( n_haps );
  for(int j = 0; j < proba_haplos.ncol(); j++) {
    for(int i = 0; i < n_haps; i++)
      p[i] = proba_haplos(i,j);
    make_inds(x, N[j], p, length_tiles); 
  }

  int n = x.size();
  List ZYG(n);
  for(int i = 0; i < n; i++) ZYG[i] = XPtr<mozza::zygote>(x[i], true);
  ZYG.attr("class") = CharacterVector::create("zygote");
  return ZYG;

}
