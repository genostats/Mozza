//[[Rcpp::depends(gaston)]]
#include <Rcpp.h>
#include "mozza.h"
#include "gaston/matrix4.h"
#include "hbd_segments.h"
using namespace Rcpp;

// cette fonction fait des individus avec des haplos 
// H et M qui partagent des segments de longueur le1 entrecoupés par des
// segments de longueur le2 (en espérance)
// Ils sont pushé dans ZYGS
// on renvoie (un vecteur de) proportion de génome HBD
std::vector<double> make_inbreds(int n, double le1, double le2, std::vector<mozza::zygote> & ZYG, int n_haps = 100, double length_tiles = 20.) {
  std::vector<double> inb;
  for(int i = 0; i < n; i++) {
    mozza::mosaic H(mozza::human_autosomes_b37, n_haps, length_tiles);
    mozza::mosaic A(mozza::human_autosomes_b37, n_haps, length_tiles);

    mozza::zygote Z( H, mozza::mosaic(H, A, le1, le2) );

    auto HBD = HBD_length( Z );
    inb.push_back( HBD / mozza::length_human_autosomes_b37  );
    ZYG.push_back(Z);
  }
  return inb;
}

//[[Rcpp::export]]
List make_inbreds(int N, double le1, double le2, double length_tiles, XPtr<matrix4> Haplos, 
                  IntegerVector chr, NumericVector dist, bool segments = false) {
  std::vector<mozza::zygote> x; 
  int n_haps = Haplos->ncol; // chaque haplotype = un "individu"
  List L;
  L["inb"] = make_inbreds(N, le1, le2, x, n_haps, length_tiles);
  mozza::mappedBed<IntegerVector, NumericVector> MB(Haplos, chr, dist);
  L["bed"] = drop_to_bed_matrix(x, MB);
  if(segments) {
    List S(N);
    for(int i = 0; i < N; i++) {
      auto seg = HBD_segments(x[i]);
      DataFrame D = DataFrame::create( Named("chr") = wrap(seg.chr), Named("beg") = wrap(seg.beg), Named("end") = wrap(seg.end) );
      S[i] = D;
    }
    L["segments"] = S;
  }
  return L;
}


