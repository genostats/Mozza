//[[Rcpp::depends(gaston)]]
#include <Rcpp.h>
#include "Mozza.h"
#include <memory>

#define _debug_mozza_forward_ false

using namespace Rcpp;

// zygotes : liste de zygotes qui constituent la première génération
// nGen : nb total de generation
// keep : nb de generations à garder
// lambda : parametre loi de Poisson (nb enfants par couples)
//[[Rcpp::export]]
List forward_(List zygotes, int nGen, int keep, double lambda) { 

  if(keep < 1) keep = 1;
  if(keep >= nGen) { // on ne peut pas garder plus de génération qu'on en a simulées...
    keep = nGen;
  }
 
  int mykeep = (keep == 1)?2:keep;  // si keep = 1 il faut garder une génération de plus pour que l'algo ne se marche pas dessus

  // vecteur de longueur keep, pour contenir les générations conservées
  std::vector<std::vector<mozza::zygote*>> POP(mykeep);

  // pour numéroter les individus et garder trace de leurs parents
  // (on aurait pu créer une classe de zygotes numérotés
  // mais il faudrait écrire tous les constructeurs...)
  std::vector<std::vector<std::tuple<int,int,int>>> NUM(mykeep);
  
  int ind = 0; // l'individu courant

  // génération 0
  for(int i = 0; i < zygotes.size(); i++) {
    // on copie les éléments de zygotes (pas envie de travailler avec des pointeurs)
    XPtr<mozza::zygote> pz = zygotes[i];
    POP[0].push_back(pz);  // Ici il y a une copie mais c'est vraiment pénible à éviter
    NUM[0].emplace_back( ++ind, 0, 0 );
  }

  // générations suivantes
  for(int gen = 1; gen < nGen; gen++) {
    int parents = (gen - 1) % mykeep;
    int enfants = gen % mykeep;
    if(gen > mykeep) { // pour ne pas appeler delete sur les éléments de 'zygotes'
      // avant de 'clear' POP, on delete les zygotes qui s'y trouvent
      if(_debug_mozza_forward_) std::cout << "index = " << enfants << ", " << POP[enfants].size() << " zygotes to delete\n";
      for(auto a : POP[enfants]) delete a;  // pour éviter les fuites de mémoire !
    }
    if(_debug_mozza_forward_) std::cout << "filling POP at index " << enfants << "\n";
    POP[enfants].clear(); 
    NUM[enfants].clear();
    int nParents = POP[parents].size();
    // on mélange les indices pour créer les couples.
    IntegerVector I = sample(nParents, nParents, false, R_NilValue, false);
    for(int k = 0; k < nParents/2; k++) {
      int nOff = R::rpois(lambda);
      for(int a = 0; a < nOff; a++) {
        mozza::zygote * pz = new mozza::zygote(*POP[parents][I[2*k]], *POP[parents][I[2*k+1]] );
        POP[enfants].push_back(pz);
        NUM[enfants].emplace_back( ++ind, std::get<0>(NUM[parents][I[2*k]]), std::get<0>(NUM[parents][I[2*k+1]]) );
      }
    }
  }

  // FINALISATION
  // On commence par récupérer les ids des individus des 'keep' dernières générations, et ceux de leurs parents
  std::vector<int> ID, FATHER, MOTHER;
  for(unsigned int gen = nGen-keep; gen < nGen; gen++) {
    int g = gen % mykeep;
    for(auto x : NUM[g]) {
      ID.push_back( std::get<0>(x) );
      FATHER.push_back( std::get<1>(x) );
      MOTHER.push_back( std::get<2>(x) );
    }
  }

  // Puis on récupère les zygotes des 'keep' dernières générations
  int N = ID.size();
  List ZYG(N);
  int k = 0;
  for(unsigned int gen = nGen-keep; gen < nGen; gen++) {
    int g = gen % mykeep;
    if(gen == 0) { // on renvoie les individus de la generation 0 (se produit si nGen == keep)
      if(_debug_mozza_forward_) std::cout << "copying XPtr from original data\n";
      // il ne faut pas refaire un XPtr vers le même objet ! (double free)
      for(int i = 0; i < zygotes.size(); i++) ZYG[k++] = zygotes[i];
    } else {
      if(_debug_mozza_forward_) std::cout << "building XPtr from individuals at index " << g << "\n";
      for(int i = 0; i < POP[g].size(); i++) ZYG[k++] = XPtr<mozza::zygote>( POP[g][i], true);
    }
    if(_debug_mozza_forward_) std::cout << "clearing POP["<< g << "] to avoid deleting zygotes sent back to user\n";
    POP[g].clear(); // il ne faudra pas 'delete' ceux là
  }

  // si la generation 0 n'a pas été overwritten, il ne faut pas la 'delete' non plus
  // ça arrive quand nGen = 2 et keep = 1 -> mykeep = 2
  // ou quand nGen = mykeep !
  if(nGen <= mykeep) {
    if(_debug_mozza_forward_) std::cout << "clearing POP[0] to avoid deleting users's zygotes\n";
    POP[0].clear();
  }
 
  // on a maintenant clear
  // et on delete les zygotes créés et non renvoyés à l'utilisateur
  // le seul cas où il y en a c'est quand keep = 1 et mykeep = 2
  // il a fallu on garder une génération de plus qu'on en renvoie
  for(unsigned int i = 0; i < mykeep; i ++) {
    if(_debug_mozza_forward_) std::cout << "index = " << i << ", " << POP[i].size() << " zygotes to delete\n";
    for(auto a : POP[i]) delete a;  
  }
 
  ZYG.attr("class") = CharacterVector::create("zygote");

  List L;

  L["zygotes"] = ZYG;
  L["N"] = N;
  L["id"] = wrap(ID);
  L["father"] = wrap(FATHER);
  L["mother"] = wrap(MOTHER);

  return L;
}
