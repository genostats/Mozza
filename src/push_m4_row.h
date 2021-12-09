#include <Rcpp.h>
using namespace Rcpp;

#ifndef PUSH_M4_ROW_TO_STD_VEC
#define PUSH_M4_ROW_TO_STD_VEC
/* essayons d'encapsuler un bout de code qui est un peu partout dans
 * Gaston... ça reservira peut-être
 */

// ncol = le nombre d'invidus / le nombre de valeurs à pusher
// cette fonction prend un vecteur d'une bed matrix et le push back dans V
template<typename T>
inline void push_m4_row(unsigned char * snp, size_t ncol, std::vector<T> & V) {
  int true_ncol = ncol/4 + ((ncol%4 == 0)?0:1);
  for(int ii = 0; ii < true_ncol-1; ii++) {
    unsigned char x = snp[ii];
    for(int ss = 0; ss < 4; ss++) {
      V.push_back(x&3);
      x >>= 2;
    }
  }
  // the last one (is there any time gain with this ?!)
  { int ii = true_ncol-1;
    unsigned char x = snp[ii];
    for(int ss = 0; ss < 4 && 4*ii+ss < ncol; ss++) {
      V.push_back(x&3);
      x >>= 2;
    }
  }
}

// et un autre bout de code : un moyen de faire du push_back sur une ligne (SNPs) de la bed matrix 
// on fait attention à ce qu'il y ait bien des NA en bout de ligne
// et on fournit un opérateur push_back qui va mettre les génotypes 1 à 1
// aucune vérification n'est faite sur le nombre de valeurs qu'on push back
class SNP_push_back {
  unsigned char * snp;
  int s;
public:
  SNP_push_back(unsigned char * snp_, int ncol) : snp(snp_), s(0) {
    int true_ncol = ncol/4 + ((ncol%4 == 0)?0:1);
    if(ncol%4 != 0) // seulement s'il y a un padding
      snp[true_ncol - 1] = (255 << ((ncol%4)*2)); // NAs en bout de lignes
  }
  inline void push_back(unsigned char a) {
    *snp &= ~(3 << s); // set to 00
    *snp |= (a << s);  // set to val
    if(s == 6) {
      s = 0;
      snp++;
    } else
      s += 2;
  }
};
#endif
