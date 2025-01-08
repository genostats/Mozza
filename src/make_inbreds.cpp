//[[Rcpp::depends(gaston)]]
#include <Rcpp.h>
#include "Mozza.h"
#include "gaston/matrix4.h"
#include "hbd_segments.h"
using namespace Rcpp;

// cette fonction fait *des pointeurs vers* des individus avec des haplos 
// H et M qui partagent des segments de longueur le1 entrecoupés par des
// segments de longueur le2 (en espérance)
// Ils sont pushé dans ZYGS
// on renvoie (un vecteur de) proportion de génome HBD
std::vector<double> make_inbreds(int n, double le1, double le2, std::vector<mozza::zygote *> & ZYG, int n_haps = 100, double length_tiles = 20.) {
  std::vector<double> inb;
  for(int i = 0; i < n; i++) {
    mozza::mosaic H(mozza::Autosomes(), n_haps, length_tiles);
    mozza::mosaic A(mozza::Autosomes(), n_haps, length_tiles);

    mozza::zygote * pZ = new mozza::zygote( H, mozza::mosaic(H, A, le1, le2) );

    auto HBD = HBD_length( *pZ );
    inb.push_back( HBD / mozza::lengthAutosomes()  );
    ZYG.push_back(pZ);
  }
  return inb;
}

// idem mais les longueurs le1 et le2 sont des vecteurs
template<typename T>
std::vector<double> make_inbreds_vecle(int n, T le1, T le2, std::vector<mozza::zygote *> & ZYG, int n_haps = 100, double length_tiles = 20.) {
  std::vector<double> inb;
  if(le1.size() != n || le2.size() != n) 
    throw std::runtime_error("Size of length tiles vectors incompatible with the number of zygotes to generate");
  for(int i = 0; i < n; i++) {
    mozza::mosaic H(mozza::Autosomes(), n_haps, length_tiles);
    mozza::mosaic A(mozza::Autosomes(), n_haps, length_tiles);

    mozza::zygote * pZ = new mozza::zygote( H, mozza::mosaic(H, A, le1[i], le2[i]) );

    auto HBD = HBD_length( *pZ );
    inb.push_back( HBD / mozza::lengthAutosomes()  );
    ZYG.push_back(pZ);
  }
  return inb;
}

//[[Rcpp::export]]
List make_inbreds(unsigned int N, unsigned int n_haps, NumericVector le1, NumericVector le2, double length_tiles, bool segments = false) {
  std::vector<mozza::zygote *> x; 
  List L;

  if(le1.size() == 1) {
    L["inbred.coef"] = make_inbreds(N, le1[0], le2[0], x, n_haps, length_tiles);
  } else {
    L["inbred.coef"] = make_inbreds_vecle(N, le1, le2, x, n_haps, length_tiles);
  }

  List ZYG(N);
  for(int i = 0; i < N; i++) ZYG[i] = XPtr<mozza::zygote>(x[i], true);
  L["zygotes"] = ZYG;

  if(segments) {
    List S(N);
    for(int i = 0; i < N; i++) {
      auto seg = HBD_segments( *x[i] );
      DataFrame D = DataFrame::create( Named("chr") = wrap(seg.chr), Named("beg") = wrap(seg.beg), Named("end") = wrap(seg.end) );
      S[i] = D;
    }
    L["segments"] = S;
  }
  return L;
}
/*
 * ye old version (which created 
 * the bed matrix)
 *
 
List make_inbreds(int N, NumericVector le1, NumericVector le2, double length_tiles, XPtr<matrix4> Haplos, 
                  IntegerVector chr, NumericVector dist, bool segments = false) {
  std::vector<mozza::zygote> x; 
  int n_haps = Haplos->ncol; // chaque haplotype = un "individu"
  List L;
  if(le1.size() == 1) {
    L["inb"] = make_inbreds(N, le1[0], le2[0], x, n_haps, length_tiles);
  } else {
    L["inb"] = make_inbreds_vecle(N, le1, le2, x, n_haps, length_tiles);
  }
  mozza::mappedBed<IntegerVector, NumericVector> MB(Haplos, chr, dist);
  L["bed"] = zygote_to_bed_matrix(x, MB);
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
*/

