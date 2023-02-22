#include <Rcpp.h>
#include <cmath>

inline double randomExp(double len) {
  if(std::isinf(len)) 
    return INFINITY;
  else
   return R::rexp(len); // R:rexp(le) tire bien dans une loi d'espérance le
}
