#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

std::string formating(double x, int pr = 3) {
  bool neg = false;
  if(x < 0) {
    neg = true;
    x = -x;
  }
  unsigned long int xx = std::round(x * std::pow(10, pr));  // si pr est trop grand on va avoir un problème
  std::string r = std::to_string(xx);
  if( r.size() <= pr ) {
    std::string a = "0.";
    for(int i = 0; i < pr - r.size(); i ++) a += "0";
    r = a + r;
  } else {
    r = r.substr(0, r.size() - pr) + "." + r.substr(r.size() - pr, pr);
  }
  if(neg) r = "-" + r;
  return r;
}

// [[Rcpp::export]]
std::string pasteVcfElts(CharacterVector genoN, IntegerVector rd, IntegerVector rd2, IntegerVector rd3, 
                         IntegerVector GQ, NumericMatrix phtab, NumericMatrix GP, NumericVector DS) {
  std::string x;
  int n = genoN.size();
  for(int i = 0; i < n; i++) {
    if(rd[i] + rd2[i] <= rd3[i]) { // "no reads"
      x += "./.:";                                                       // GT
      x += std::to_string(rd[i]) + "," + std::to_string(rd2[i]) + ":";   // AD
      x += std::to_string(rd[i]+rd2[i]+rd3[i]) + ":";                    // DP
      x += ".:";                                                         // GQ
      x += ".,.,.:";                                                     // PL
      x += ".,.,.:";                                                     // GP
      x += ".";                                                          // DS
    } else {
      x += std::string(CHAR(genoN[i])) + ":";                            // GT
      x += std::to_string(rd[i]) + "," + std::to_string(rd2[i]) + ":";   // AD
      x += std::to_string(rd[i]+rd2[i]+rd3[i]) + ":";                    // DP
      x += std::to_string(GQ[i]) + ":";                                  // GQ
      x += std::to_string((int) round(phtab(i,0))) + ",";                // PL
      x += std::to_string((int) round(phtab(i,1))) + ",";
      x += std::to_string((int) round(phtab(i,2))) + ":"; 
      x += formating(GP(i,0)) + "," + formating(GP(i,1)) + "," ;         // GP
      x += formating(GP(i,2)) + ":";
      x += formating(DS(i));                                             // DS
    }
    if(i != n-1) x += "\t";
  } 
  return x;
}
