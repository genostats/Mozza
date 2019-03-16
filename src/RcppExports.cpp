// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/Mozza.h"
#include <Rcpp.h>

using namespace Rcpp;

// MH_cpp
NumericMatrix MH_cpp(int B, double sd, int burn, int thin);
RcppExport SEXP _Mozza_MH_cpp(SEXP BSEXP, SEXP sdSEXP, SEXP burnSEXP, SEXP thinSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type B(BSEXP);
    Rcpp::traits::input_parameter< double >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    rcpp_result_gen = Rcpp::wrap(MH_cpp(B, sd, burn, thin));
    return rcpp_result_gen;
END_RCPP
}
// make_inds
List make_inds(int n, double length_tiles, XPtr<matrix4> Haplos, IntegerVector chr, NumericVector dist, bool kinship);
RcppExport SEXP _Mozza_make_inds(SEXP nSEXP, SEXP length_tilesSEXP, SEXP HaplosSEXP, SEXP chrSEXP, SEXP distSEXP, SEXP kinshipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type length_tiles(length_tilesSEXP);
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type Haplos(HaplosSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type chr(chrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dist(distSEXP);
    Rcpp::traits::input_parameter< bool >::type kinship(kinshipSEXP);
    rcpp_result_gen = Rcpp::wrap(make_inds(n, length_tiles, Haplos, chr, dist, kinship));
    return rcpp_result_gen;
END_RCPP
}
// make_pairs
List make_pairs(int N, double le1, double le2, double length_tiles, XPtr<matrix4> Haplos, IntegerVector chr, NumericVector dist);
RcppExport SEXP _Mozza_make_pairs(SEXP NSEXP, SEXP le1SEXP, SEXP le2SEXP, SEXP length_tilesSEXP, SEXP HaplosSEXP, SEXP chrSEXP, SEXP distSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type le1(le1SEXP);
    Rcpp::traits::input_parameter< double >::type le2(le2SEXP);
    Rcpp::traits::input_parameter< double >::type length_tiles(length_tilesSEXP);
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type Haplos(HaplosSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type chr(chrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dist(distSEXP);
    rcpp_result_gen = Rcpp::wrap(make_pairs(N, le1, le2, length_tiles, Haplos, chr, dist));
    return rcpp_result_gen;
END_RCPP
}
// nuclear_families
List nuclear_families(int nb_fams, int nb_offsprings, double tile_length, XPtr<matrix4> Haplos, IntegerVector chr, NumericVector dist, bool kinship);
RcppExport SEXP _Mozza_nuclear_families(SEXP nb_famsSEXP, SEXP nb_offspringsSEXP, SEXP tile_lengthSEXP, SEXP HaplosSEXP, SEXP chrSEXP, SEXP distSEXP, SEXP kinshipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nb_fams(nb_famsSEXP);
    Rcpp::traits::input_parameter< int >::type nb_offsprings(nb_offspringsSEXP);
    Rcpp::traits::input_parameter< double >::type tile_length(tile_lengthSEXP);
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type Haplos(HaplosSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type chr(chrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dist(distSEXP);
    Rcpp::traits::input_parameter< bool >::type kinship(kinshipSEXP);
    rcpp_result_gen = Rcpp::wrap(nuclear_families(nb_fams, nb_offsprings, tile_length, Haplos, chr, dist, kinship));
    return rcpp_result_gen;
END_RCPP
}
// essai
void essai();
RcppExport SEXP _Mozza_essai() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    essai();
    return R_NilValue;
END_RCPP
}
// test_cursor
void test_cursor();
RcppExport SEXP _Mozza_test_cursor() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    test_cursor();
    return R_NilValue;
END_RCPP
}
// essai2
void essai2();
RcppExport SEXP _Mozza_essai2() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    essai2();
    return R_NilValue;
END_RCPP
}
// essai3
void essai3();
RcppExport SEXP _Mozza_essai3() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    essai3();
    return R_NilValue;
END_RCPP
}
// essai3bis
void essai3bis();
RcppExport SEXP _Mozza_essai3bis() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    essai3bis();
    return R_NilValue;
END_RCPP
}
// test_sharing
List test_sharing(int n, double le1, double le2);
RcppExport SEXP _Mozza_test_sharing(SEXP nSEXP, SEXP le1SEXP, SEXP le2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type le1(le1SEXP);
    Rcpp::traits::input_parameter< double >::type le2(le2SEXP);
    rcpp_result_gen = Rcpp::wrap(test_sharing(n, le1, le2));
    return rcpp_result_gen;
END_RCPP
}
// test_IBD_sibs
List test_IBD_sibs(int n, int n_haps, double length_tiles);
RcppExport SEXP _Mozza_test_IBD_sibs(SEXP nSEXP, SEXP n_hapsSEXP, SEXP length_tilesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type n_haps(n_hapsSEXP);
    Rcpp::traits::input_parameter< double >::type length_tiles(length_tilesSEXP);
    rcpp_result_gen = Rcpp::wrap(test_IBD_sibs(n, n_haps, length_tiles));
    return rcpp_result_gen;
END_RCPP
}
// test_push_genos
std::vector<int> test_push_genos(IntegerVector H);
RcppExport SEXP _Mozza_test_push_genos(SEXP HSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type H(HSEXP);
    rcpp_result_gen = Rcpp::wrap(test_push_genos(H));
    return rcpp_result_gen;
END_RCPP
}
// test_xptr
XPtr<matrix4> test_xptr(XPtr<matrix4> Haplos, IntegerVector chr, NumericVector dist);
RcppExport SEXP _Mozza_test_xptr(SEXP HaplosSEXP, SEXP chrSEXP, SEXP distSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type Haplos(HaplosSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type chr(chrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dist(distSEXP);
    rcpp_result_gen = Rcpp::wrap(test_xptr(Haplos, chr, dist));
    return rcpp_result_gen;
END_RCPP
}
// families_of_4_v0
XPtr<matrix4> families_of_4_v0(int N, XPtr<matrix4> Haplos, IntegerVector chr, NumericVector dist);
RcppExport SEXP _Mozza_families_of_4_v0(SEXP NSEXP, SEXP HaplosSEXP, SEXP chrSEXP, SEXP distSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type Haplos(HaplosSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type chr(chrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dist(distSEXP);
    rcpp_result_gen = Rcpp::wrap(families_of_4_v0(N, Haplos, chr, dist));
    return rcpp_result_gen;
END_RCPP
}
// families_of_4
XPtr<matrix4> families_of_4(int N, XPtr<matrix4> Haplos, IntegerVector chr, NumericVector dist);
RcppExport SEXP _Mozza_families_of_4(SEXP NSEXP, SEXP HaplosSEXP, SEXP chrSEXP, SEXP distSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type Haplos(HaplosSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type chr(chrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dist(distSEXP);
    rcpp_result_gen = Rcpp::wrap(families_of_4(N, Haplos, chr, dist));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Mozza_MH_cpp", (DL_FUNC) &_Mozza_MH_cpp, 4},
    {"_Mozza_make_inds", (DL_FUNC) &_Mozza_make_inds, 6},
    {"_Mozza_make_pairs", (DL_FUNC) &_Mozza_make_pairs, 7},
    {"_Mozza_nuclear_families", (DL_FUNC) &_Mozza_nuclear_families, 7},
    {"_Mozza_essai", (DL_FUNC) &_Mozza_essai, 0},
    {"_Mozza_test_cursor", (DL_FUNC) &_Mozza_test_cursor, 0},
    {"_Mozza_essai2", (DL_FUNC) &_Mozza_essai2, 0},
    {"_Mozza_essai3", (DL_FUNC) &_Mozza_essai3, 0},
    {"_Mozza_essai3bis", (DL_FUNC) &_Mozza_essai3bis, 0},
    {"_Mozza_test_sharing", (DL_FUNC) &_Mozza_test_sharing, 3},
    {"_Mozza_test_IBD_sibs", (DL_FUNC) &_Mozza_test_IBD_sibs, 3},
    {"_Mozza_test_push_genos", (DL_FUNC) &_Mozza_test_push_genos, 1},
    {"_Mozza_test_xptr", (DL_FUNC) &_Mozza_test_xptr, 3},
    {"_Mozza_families_of_4_v0", (DL_FUNC) &_Mozza_families_of_4_v0, 4},
    {"_Mozza_families_of_4", (DL_FUNC) &_Mozza_families_of_4, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_Mozza(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
