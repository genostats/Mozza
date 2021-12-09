#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include "gaston/matrix4.h"
#include "push_m4_row.h"

// [[Rcpp::export]]
void write_hap_file(XPtr<matrix4> p_A, std::string filename) {
  std::ofstream file(filename, std::ofstream::binary);
  if(!file.is_open()) {
    stop("Cannot open file");
  }

  int n_snp = p_A->nrow;
  int n_ind = p_A->ncol;
  std::vector<char> v;
  char sp = 32;
  char lf = 10;
  for(int i = 0; i < n_snp; i++) {
    // on recupere un vecteur de génotypes
    v.clear();
    v.reserve(n_ind);
    push_m4_row(p_A->data[i], n_ind, v);
    for(int j = 0; j < n_ind; j++) {
      char x = '0' + v[j];
      file.write(&x, 1);
      if(j + 1 == n_ind) 
        file.write(&lf, 1);
      else
        file.write(&sp, 1);
    }
  }

  file.close();
}

