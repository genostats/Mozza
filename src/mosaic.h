#ifndef MOSAIC_
#define MOSAIC_
#define SHOW(x) Rcpp::Rcout << #x << " = " << (x) << std::endl;

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
  mosaic(const std::vector<double> & chr_len, int tile = 0);
  
  // constructeur par mélange 
  // le1, le2 = longueur moyenne des morceaux (loi exp) des deux haplotypes d'origine
  // sauf erreur les valeurs par défaut correspondent à la meiose (morceaux d'un morgan)
  mosaic(mosaic & M1, mosaic & M2, double le1 = 100., double l2 = 100.);
  
  // constructeur avec des tuiles tirée uniforméments entre 0 et (ntiles - 1), 
  // de longueurs prises dans une loi exp lambda = 1/mean_length
  // [truc à la hapgen]
  mosaic(const std::vector<double> & chr_len, int ntiles, double mean_length);
  
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

};
#endif