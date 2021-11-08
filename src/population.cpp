//[[Rcpp::depends(gaston)]]
#include <Rcpp.h>
#include "mozza.h"
using namespace Rcpp;

// n0 : nb indiv à la génération 0
// nGen : nb total de generation
// keep : nb de generations à garder
// lambda : parametre loi de Poisson
//[[Rcpp::export]]
List population(int n0, int nGen, int keep, double lambda, 
                double tile_length, XPtr<matrix4> Haplos, IntegerVector chr, NumericVector dist, 
                bool kinship = false, bool fraternity = false) {

  if(keep < 1) keep = 1;
  // vecteur de longueur keep, pour contenir les générations conservées
  std::vector<std::vector<mozza::zygote>> POP(keep);
  // et pour numéroter les individus et garder trace de leurs parents
  // (on aurait pu créer une classe de zygotes numérotés
  // mais il faudrait écrire tous les constructeurs...)
  std::vector<std::vector<std::tuple<int,int,int>>> NUM(keep);
  
  int n_haps = Haplos->ncol; // chaque haplotype = un "individu" de la bed matrix
  int ind = 0; // l'individu courant
  // génération 0
  for(int i = 0; i < n0; i++) {
    POP[0].emplace_back( mozza::human_autosomes_b37, n_haps, tile_length ); // appelle le constructeur de mozza::zygote
    NUM[0].emplace_back( ++ind, 0, 0 );
  }

  // générations suivantes
  for(int gen = 1; gen < nGen; gen++) {
    int parents = (gen - 1) % keep;
    int enfants = gen % keep;
    POP[enfants].clear(); 
    NUM[enfants].clear();
    int nParents = POP[parents].size();
    // on mélange les indices pour créer les couples.
    IntegerVector I = sample(nParents, nParents, false, R_NilValue, false);
    for(int k = 0; k < nParents/2; k++) {
      int nOff = R::rpois(lambda);
      for(int a = 0; a < nOff; a++) {
        POP[enfants].push_back( POP[parents][I[2*k]] + POP[parents][I[2*k+1]] );
        NUM[enfants].emplace_back( ++ind, std::get<0>(NUM[parents][I[2*k]]), std::get<0>(NUM[parents][I[2*k+1]]) );
      }
    }
  }

  // Une fois que c'est fini on met tout dans un seul vecteur.
  std::vector<mozza::zygote> ZYG;
  // Et on récupère les ids, et ceux des parents
  std::vector<int> ID, FATHER, MOTHER;
  for(int gen = nGen-keep; gen < nGen; gen++) {
    int g = gen % keep;
    for(auto zy : POP[g]) {
      ZYG.push_back(zy);
    }
    for(auto x : NUM[g]) {
      ID.push_back( std::get<0>(x) );
      FATHER.push_back( std::get<1>(x) );
      MOTHER.push_back( std::get<2>(x) );
    }
  }
  List L;

  L["bed"] = drop_to_bed_matrix(ZYG, Haplos, chr, dist);
  L["N"] = ZYG.size();
  L["id"] = wrap(ID);
  L["father"] = wrap(FATHER);
  L["mother"] = wrap(MOTHER);

  if(kinship) 
    L["kinship"] = kinship_matrix(ZYG);
  if(fraternity) 
    L["fraternity"] = fraternity_matrix(ZYG);

  return L;
}
/*
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

// ce code ressemble à celui de make_inds_probs
//
// proba_haplos est une matrice, chaque colonnes donne un jeu de proba sur les haplotypes
// le vecteur d'effectif N contient les effectifs à générer pour chacune des colonnees
// Cela permet de générer plus facilement des données avec des sous-populations ayant des 
// proportions différentes.
//
// Nfams = autant d'élts que de colonnes à proba_haplos
// chaque famille a nb_offsprings
//[[Rcpp::export]]
List nuclear_families_probs(IntegerVector Nfams, int nb_offsprings, NumericMatrix proba_haplos, double length_tiles,
                     XPtr<matrix4> Haplos, IntegerVector chr, NumericVector dist, bool kinship = false, bool fraternity = false) {

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
  if(fraternity)
    L["fraternity"] = fraternity_matrix(ZYG);
  return L;
}
*/
