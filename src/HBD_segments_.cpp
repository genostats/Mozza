#include <Rcpp.h>
#include "Mozza.h"
#include "hbd_segments.h"

//[[Rcpp::export]]
Rcpp::List HBD_segments_(Rcpp::List L) {
  int N = L.size();
  Rcpp::List S(N);
  for(int i = 0; i < N; i++) {
    Rcpp::XPtr<mozza::zygote> pz = L[i];
    auto seg = HBD_segments(*pz, 1);
    Rcpp::DataFrame D = Rcpp::DataFrame::create( Named("chr") = wrap(seg.chr), Named("beg") = wrap(seg.beg), Named("end") = wrap(seg.end) );
    S[i] = D;
  }
  return S;
}

