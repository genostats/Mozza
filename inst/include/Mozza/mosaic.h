#ifndef MOSAIC_
#define MOSAIC_
#define SHOW(x) Rcpp::Rcout << #x << " = " << (x) << std::endl;

#include "randomExp.h"
#include "sampleInt.h"
#include "sample_cp.h"
/*
 * les tuiles de la mosaique sont fermées à gauche [ --- ] ( --- ]  ...  ( --- ]
 * sauf la première (la position 0 est dedans) 
 * 
 * La dernière tuile d'un chr doit s'arrêter sur chr_len[chr]
 * 
 * Les vecteurs de bpoints donnent les bornes à droite des tuiles 
 * (le début en 0 est implicite, la dernière tuile d'un chromomome 
 * finit sur chr_len[chr]) 
 * 
 * Dans le cas général quand cursor_pos est pile sur la position d'un bpoint, 
 * i_cursor pointe "la tuile qui est à gauche" (et donc bpoint[i_cursor] == cursor_pos) 
 * puisqu'on considère que ce point est dans la tuile de gauche. 
 * 
 * Cf fonction test_cursor pour un exemple. 
 */
namespace mozza {
class mosaic {
public:
  int chrs; // nombre de chromosomes (eg 22)
  // référence à un vecteur de chrs (22) longueurs (en cM...)
  const std::vector<double> & chr_len; 
  
  // longueur totale des chromosomes
  double genome_length;
  
  // un vecteur de chrs (22) vecteurs de 'tiles'
  std::vector<std::vector<int>> tiles; 

  // idem avec les points de coupure (en cM) entre les tuiles
  // NB Il y a toujours un bpoint en fin de chr, sa position doit coincider avec chr_len[chr]
  std::vector<std::vector<double>> bpoints; 

  // position du curseur
  unsigned int cursor_chr;  
  double cursor_pos; // en cM
  unsigned int i_cursor;
  
  // constructeur avec la tuile tile (0) répétée sur tous les chrs
  mosaic(const std::vector<double> & chr_len, int tile = 0) : 
      chrs(chr_len.size()), chr_len(chr_len), 
      genome_length(std::accumulate(chr_len.begin(), chr_len.end(), 0.)),
      tiles(chrs), bpoints(chrs), cursor_chr(0), cursor_pos(0.), i_cursor(0) {
    for(int i = 0; i < chrs; i++) {
      tiles[i].push_back(tile);
      bpoints[i].push_back(chr_len[i]);
    }
  }

  
  // constructeur avec des tuiles tirée uniforméments entre 0 et (ntiles - 1), 
  // de longueurs prises dans une loi exp lambda = 1/mean_length
  // [truc à la hapgen]
  mosaic(const std::vector<double> & chr_len, int ntiles, double mean_length) :
      chrs(chr_len.size()), chr_len(chr_len), 
      genome_length(std::accumulate(chr_len.begin(), chr_len.end(), 0.)),
      tiles(chrs), bpoints(chrs), cursor_chr(0), cursor_pos(0.), i_cursor(0) {
    for(int i = 0; i < chrs; i++) {
      double le = chr_len[i];
      double pos = 0;
      while(true) {
        // tire et ajoute une tuile
        tiles[i].push_back(R::runif(0,1)*ntiles);
        // position de fin de cette tuile
        pos += randomExp(mean_length);
        if(pos < le) {
          bpoints[i].push_back(pos);
        } else {
          bpoints[i].push_back(le);
          break;
        }
      }
    }
  }

 
  // constructeur avec des tuiles de 0 à proba_tiles.size()-1, 
  // de longueur prise dans une exponentielle comme ci-dessus
  // d'indices tirés selon les valeurs du vecteur proba_tiles
  template<typename vec> 
  mosaic(const std::vector<double> & chr_len, const vec & proba_tiles, double mean_length) :
      chrs(chr_len.size()), chr_len(chr_len), genome_length(std::accumulate(chr_len.begin(), chr_len.end(), 0.)),
      tiles(chrs), bpoints(chrs), cursor_chr(0), cursor_pos(0.), i_cursor(0) {
    // probas cumulées
    std::vector<double> cum_probas;
    double S = 0;
    for(auto a : proba_tiles) {
      S += a;
      cum_probas.push_back(S);
    }
    // normaliser !
    for(auto & a : cum_probas) a /= S;
    // code quasi identique à celui du constructeur précédent
    for(int i = 0; i < chrs; i++) {
      double le = chr_len[i];
      double pos = 0;
      while(true) {
        // tire et ajoute une tuile
          tiles[i].push_back( sample_cp(cum_probas) );
        // position de fin de cette tuile
        pos += randomExp(mean_length);
        if(pos < le) {
          bpoints[i].push_back(pos);
        } else {
          bpoints[i].push_back(le);
          break;
        }
      }
    }
  }
  
  // constructeur par mélange 
  // le1, le2 = longueur moyenne des morceaux (loi exp) des deux haplotypes d'origine
  // sauf erreur les valeurs par défaut correspondent à la meiose (morceaux d'un morgan)
  mosaic(mosaic & M1, mosaic & M2, double le1 = 100., double le2 = 100.) :
      chrs(M1.chrs), chr_len(M1.chr_len), genome_length(M1.genome_length),
      tiles(chrs),  bpoints(chrs), cursor_chr(0), cursor_pos(0.), i_cursor(0) {

    double p1 = le1/(le1 + le2); // proba d'être sur M1.
    if(std::isnan(p1)) {
      if(std::isinf(le1) && !std::isinf(le2))
        p1 = double(1);
      else
        throw std::invalid_argument("nan lengthes or two infinite lengthes");
    }
    if(M2.chrs != chrs)
      throw std::invalid_argument("Two mosaic haplotypes with different chrs");
    for(int i = 0; i < chrs; i++) {
      M1.set_cursor(i);
      M2.set_cursor(i);

      if(chr_len[i] != M2.chr_len[i])
        throw std::invalid_argument("Chromosomes of different lengths");
      // on tire sur quel haplo on commence
      bool on_m1 = (R::runif(0,1) < p1);
      // on avance...
      while(M1.cursor_pos < chr_len[i]) {
        double le = on_m1?randomExp(le1):randomExp(le2);
        if(on_m1) {
          M1.step_cursor(le, tiles[i], bpoints[i]);
          M2.step_cursor(le);
        } else {
          M1.step_cursor(le);
          M2.step_cursor(le, tiles[i], bpoints[i]);
        }
        // il faut ajouter le point de crossover
        if(on_m1)
          tiles[i].push_back(M1.tile_at_cursor());
        else
          tiles[i].push_back(M2.tile_at_cursor());
        bpoints[i].push_back(M1.cursor_pos);
  
        on_m1 = !on_m1; // on change
      }
    }
  }

  // positionner le curseur
  void set_cursor(unsigned int chr = 0, double pos = 0.);

  // avancer le curseur [en récupérant ou pas les tuiles et les break points rencontrées sur le chemin]
  void step_cursor(double length); // avance de length
  void forward_cursor(double pos); // avance jusqu'à pos (supposé > cursor_pos, sinon on bouge pas)

  // idem step cursor mais push_back dans til et pbt les numéros de tuiles et les positions de bpoints
  // "dépassés sur le chemin"
  void step_cursor(double length, std::vector<int> & til, std::vector<double> & bpt);

  // numéro de la tuile au curseur
  int tile_at_cursor() {
    return tiles[cursor_chr][i_cursor];
  }
  
  // à des fins de debugage...
  void print();
  void print_chr(int i);
  void print_cursor();
};

}
#endif
