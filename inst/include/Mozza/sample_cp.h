#include <Rcpp.h>
#ifndef __sample_cp__
#define __sample_cp__

namespace mozza {

// on donne un vecteur de proba cumulées
// exemple : 0.1, 0.4, 1
// -> renvoie 0 avec proba 0.1, 2 avec proba 0.3, 4 avec proba 0.6
inline int sample_cp(const std::vector<double> & cum_prob) {
  double r = unif_rand();
  int i = 0;
  while(true) {
    if(r <= cum_prob[i])
      break;
    ++i;
  }
  return i;
}

} // namespace mozza

#endif

