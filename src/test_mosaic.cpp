//[[Rcpp::depends(gaston)]]
#include <Rcpp.h>
#include "Mozza.h"
#include "gaston/matrix4.h"
using namespace Rcpp;

// [[Rcpp::export]]
void essai() {
  mozza::mosaic M(mozza::human_autosomes_b37, 10, 50.0);
  // Rcout << M.chr_len[0] << "\n";
  // Rcout << M.bpoints[0][0] << "\n";
  // Rcout << M.tiles[0][0] << "\n";
  // Rcout << M.tiles.size() << "\n";
  mozza::mosaic L(M);
  L.print();
  M.set_cursor(9);
  std::vector<int> til;
  std::vector<double> pos;
  for(int i = 0; i < 10 ; i++) {
    M.step_cursor(10);
    Rcout << M.cursor_pos << ":" << M.tile_at_cursor() << " ";
    /* Rcout << M.cursor_pos << ":\n";
    M.step_cursor(10, til, pos);
    for(auto & x : til) Rcout << x << " "; Rcout << "\n";
    for(auto & x : pos) Rcout << x << " "; Rcout << "\n";
    */
  }
  Rcout << "\n";
}

//[[Rcpp::export]]
void test_cursor() {
  mozza::mosaic M(mozza::human_autosomes_b37, 10, 50.0);
  M.print_chr(0);
  Rcpp::Rcout << "set_cursor(0)\n";
  M.set_cursor(0);
  SHOW( M.cursor_chr );
  SHOW( M.cursor_pos );
  SHOW( M.i_cursor );
  SHOW( M.bpoints[M.cursor_chr][M.i_cursor] );
  SHOW( M.tile_at_cursor() );

  Rcpp::Rcout << "\nstep_cursor to bpoint\n";
  M.step_cursor(  M.bpoints[M.cursor_chr][M.i_cursor] );
  SHOW( M.cursor_chr );
  SHOW( M.cursor_pos );
  SHOW( M.i_cursor );
  SHOW( M.bpoints[M.cursor_chr][M.i_cursor] );
  SHOW( M.tile_at_cursor() );

  Rcpp::Rcout << "\nset_cursor to bpoint\n";
  M.set_cursor(0,  M.bpoints[M.cursor_chr][M.i_cursor] );
  SHOW( M.cursor_chr );
  SHOW( M.cursor_pos );
  SHOW( M.i_cursor );
  SHOW( M.bpoints[M.cursor_chr][M.i_cursor] );
  SHOW( M.tile_at_cursor() );
  
  Rcpp::Rcout << "\nstep_cursor by 0.0001\n";
  M.step_cursor( 0.0001 );
  SHOW( M.cursor_chr );
  SHOW( M.cursor_pos );
  SHOW( M.i_cursor );
  SHOW( M.bpoints[M.cursor_chr][M.i_cursor] );
  SHOW( M.tile_at_cursor() );
  
  Rcpp::Rcout << "\nset_cursor to (0,chr_len[0])\n";
  M.set_cursor(0,  M.chr_len[0] );
  SHOW( M.cursor_chr );
  SHOW( M.cursor_pos );
  SHOW( M.i_cursor );
  SHOW( M.bpoints[M.cursor_chr][M.i_cursor] );
  SHOW( M.tile_at_cursor() );
  
  Rcpp::Rcout << "\nstep_cursor by 0.0001\n";
  M.step_cursor( 0.0001 );
  SHOW( M.cursor_chr );
  SHOW( M.cursor_pos );
  SHOW( M.i_cursor );
  SHOW( M.bpoints[M.cursor_chr][M.i_cursor] );
  SHOW( M.tile_at_cursor() );
  
}

// [[Rcpp::export]]
void essai2() {
  mozza::mosaic M1(mozza::human_autosomes_b37, 10, 20);
  mozza::mosaic M2(mozza::human_autosomes_b37, 10, 20);
  mozza::mosaic M(M1, M2);
  M1.print_chr(1);
  M2.print_chr(1);
  M.print_chr(1);
}

// [[Rcpp::export]]
void essai3() {
  mozza::mosaic M1(mozza::human_autosomes_b37, 1);
  mozza::mosaic M2(mozza::human_autosomes_b37, 2);
  mozza::mosaic M(M1, M2, 100, 100);
  Rcpp::Rcout << "M  et M1 " << mozza::IBD_sharing(M1, M) << "\n";
  Rcpp::Rcout << "M  et M2 " << mozza::IBD_sharing(M2, M) << "\n";
  // Rcpp::Rcout << "M1 et M2 " << mozza::IBD_sharing(M1, M2) << "\n";
}

// [[Rcpp::export]]
void essai3bis() {
  std::vector<double> len {293.4}; // un seul chr (idem chr 1)
  mozza::mosaic M1(len, 1);
  mozza::mosaic M2(len, 2);; 
  mozza::mosaic M(M1, M2, 100, 100);
  M1.print_chr(0);
  M2.print_chr(0);
  M.print_chr(0);
  Rcpp::Rcout << "M  et M1 " << mozza::IBD_sharing(M1, M) << "\n";
  Rcpp::Rcout << "M  et M2 " << mozza::IBD_sharing(M2, M) << "\n";
}

// [[Rcpp::export]]
List test_IBD_sharing(int n, double le1 = 100., double le2 = 100.) {
  std::vector<double> R1, R2, R3;
  for(int i = 0; i < n; i++) {
    mozza::mosaic M1(mozza::human_autosomes_b37, 100, 20);
    mozza::mosaic M2(mozza::human_autosomes_b37, 100, 20);
    // mozza::mosaic M1(mozza::human_autosomes_b37, 1);
    // mozza::mosaic M2(mozza::human_autosomes_b37, 2);
    mozza::mosaic M(M1, M2, le1, le2);
    R1.push_back( mozza::IBD_sharing(M,M1)  / mozza::length_human_autosomes_b37);
    R2.push_back( mozza::IBD_sharing(M,M2)  / mozza::length_human_autosomes_b37);
    R3.push_back( mozza::IBD_sharing(M1,M2) / mozza::length_human_autosomes_b37);    
  }
  List L;
  L["S1"] = R1;
  L["S2"] = R2;
  L["S12"] = R3;
  return L;
}

// [[Rcpp::export]]
List test_IBD_sibs(int n, int n_haps = 100, double length_tiles = 20.) {
  std::vector<double> R0, R1, R2;
  for(int i = 0; i < n; i++) {
    mozza::mosaic M1(mozza::human_autosomes_b37, n_haps, length_tiles);
    mozza::mosaic M2(mozza::human_autosomes_b37, n_haps, length_tiles);
    mozza::mosaic P1(mozza::human_autosomes_b37, n_haps, length_tiles);
    mozza::mosaic P2(mozza::human_autosomes_b37, n_haps, length_tiles);
    mozza::zygote Z1( mozza::mosaic(P1, P2), mozza::mosaic(M1, M2) );
    mozza::zygote Z2( mozza::mosaic(P1, P2), mozza::mosaic(M1, M2) );

    auto IBD = IBD_length( Z1, Z2 );
    
    R0.push_back( std::get<0>(IBD) / mozza::length_human_autosomes_b37);
    R1.push_back( std::get<1>(IBD) / mozza::length_human_autosomes_b37);
    R2.push_back( std::get<2>(IBD) / mozza::length_human_autosomes_b37);
  }
  List L;
  L["IBD0"] = R0;
  L["IBD1"] = R1;
  L["IBD2"] = R2;
  return L;
}

// H de longueur 4 au moins !! c(1L, 10L, 100L, 1000L) permet de bien vérifier
//[[Rcpp::export]]
std::vector<int> test_push_genos(IntegerVector H) {
  mozza::mosaic M1(mozza::human_autosomes_b37, 0);
  mozza::mosaic M2(mozza::human_autosomes_b37, 1);
  mozza::mosaic P1(mozza::human_autosomes_b37, 2);
  mozza::mosaic P2(mozza::human_autosomes_b37, 3);
  mozza::mosaic Ma1(P1, P2, 100, 100);
  mozza::mosaic Ma2(M1, M2, 100, 100);
  mozza::mosaic Mb1(P1, P2, 100, 100);
  mozza::mosaic Mb2(M1, M2, 100, 100);
  Ma1.print_chr(0);
  Ma2.print_chr(0);
  Mb1.print_chr(0);
  Mb2.print_chr(0);
  std::vector<mozza::zygote> x 
  { mozza::zygote(M1,  M2),
    mozza::zygote(P1,  P2),
    mozza::zygote(Ma1, Ma2), 
    mozza::zygote(Mb1, Mb2) };
  std::vector<int> R;
  push_genotypes_at_cursor(x, H, R);
  
  return R;
}

//[[Rcpp::export]]
XPtr<matrix4> test_xptr(XPtr<matrix4> Haplos, IntegerVector chr, NumericVector dist) {
  mozza::mosaic M1(mozza::human_autosomes_b37, 0);
  mozza::mosaic M2(mozza::human_autosomes_b37, 1);
  mozza::mosaic P1(mozza::human_autosomes_b37, 2);
  mozza::mosaic P2(mozza::human_autosomes_b37, 3);
  mozza::mosaic Ma1(P1, P2, 100, 100);
  mozza::mosaic Ma2(M1, M2, 100, 100);
  mozza::mosaic Mb1(P1, P2, 100, 100);
  mozza::mosaic Mb2(M1, M2, 100, 100);
  std::vector<mozza::zygote> x { 
    mozza::zygote(M1,  M2),  mozza::zygote(P1,  P2),
    mozza::zygote(Ma1, Ma2), mozza::zygote(Mb1, Mb2) };
 
  mozza::mappedBed<IntegerVector, NumericVector> MB(Haplos, chr, dist);
  return zygote_to_bed_matrix(x, MB);
}

//[[Rcpp::export]]
XPtr<matrix4> families_of_4_v0(int N, XPtr<matrix4> Haplos, IntegerVector chr, NumericVector dist) {
  std::vector<mozza::zygote> x; 
  int n_haps = Haplos->ncol; // chaque haplotype = un "individu"

  for(int i = 0; i < N; i++) {
    mozza::mosaic M1(mozza::human_autosomes_b37, n_haps, 20);
    mozza::mosaic M2(mozza::human_autosomes_b37, n_haps, 20);
    mozza::mosaic P1(mozza::human_autosomes_b37, n_haps, 20);
    mozza::mosaic P2(mozza::human_autosomes_b37, n_haps, 20);
    mozza::mosaic Ma1(P1, P2, 100, 100);
    mozza::mosaic Ma2(M1, M2, 100, 100);
    mozza::mosaic Mb1(P1, P2, 100, 100);
    mozza::mosaic Mb2(M1, M2, 100, 100);
    x.push_back(mozza::zygote(M1,  M2));
    x.push_back(mozza::zygote(P1,  P2));
    x.push_back(mozza::zygote(Ma1, Ma2)); 
    x.push_back(mozza::zygote(Mb1, Mb2));
    /* if(i == 0) {
      SHOW(mozza::IBD_sharing(M1, M2) / mozza::length_human_autosomes_b37);
      SHOW(mozza::IBD_sharing(P1, P2) / mozza::length_human_autosomes_b37);
      SHOW(mozza::IBD_sharing(P1, Ma1) / mozza::length_human_autosomes_b37);
      SHOW(mozza::IBD_sharing(Ma1,Mb1&) / mozza::length_human_autosomes_b37);
    }*/
  }
  mozza::mappedBed<IntegerVector, NumericVector> MB(Haplos, chr, dist);
  return zygote_to_bed_matrix(x, MB);
}


// new version using new constructor / + operator
//[[Rcpp::export]]
XPtr<matrix4> families_of_4(int N, XPtr<matrix4> Haplos, IntegerVector chr, NumericVector dist) {
  std::vector<mozza::zygote> x; 
  int n_haps = Haplos->ncol; // chaque haplotype = un "individu"
  
  for(int i = 0; i < N; i++) {
    mozza::zygote M(mozza::human_autosomes_b37, n_haps, 20);
    mozza::zygote F(mozza::human_autosomes_b37, n_haps, 20);
    x.push_back(M);
    x.push_back(F);
    x.push_back(M+F);
    x.push_back(M+F);
  }
  mozza::mappedBed<IntegerVector, NumericVector> MB(Haplos, chr, dist);
  return zygote_to_bed_matrix(x, MB);
}
