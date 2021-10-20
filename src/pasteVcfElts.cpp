#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::string pasteVcfElts(CharacterVector genoN, IntegerVector rd, IntegerVector rd2, IntegerVector rd3, 
                         IntegerVector GQ, NumericMatrix phtab) {
  std::string x;
  int n = genoN.size();
  for(int i = 0; i < n; i++) {
    if(rd[i] + rd2[i] <= rd3[i]) { // "no reads"
      x += "./.:";
      x += std::to_string(rd[i]) + "," + std::to_string(rd2[i]) + ":";
      x += std::to_string(rd[i]+rd2[i]+rd3[i]) + ":";
      x += ".:";
      x += ".,.,.";
    } else {
      x += CHAR(genoN[i]);
      x += ":"; 
      x += std::to_string(rd[i]) + "," + std::to_string(rd2[i]) + ":";
      x += std::to_string(rd[i]+rd2[i]+rd3[i]) + ":";
      x += std::to_string(GQ[i]) + ":";
      x += std::to_string((int) round(phtab(i,0))) + ","; 
      x += std::to_string((int) round(phtab(i,1))) + ",";
      x += std::to_string((int) round(phtab(i,2))); 
    }
    if(i != n-1) x += "\t";
  } 
  return x;
}
