#include <Rcpp.h>
#include "metropolis.h"
#include "Pi.h"

using namespace Rcpp;

//[[Rcpp::export]]
NumericMatrix MH_cpp(int B, double sd, int burn = 0, int thin = 1) {
  Pi_family F(1.0, 0.5, 1.0, -3.0);
  return Metropolis(F, B, sd, burn, thin);
}
