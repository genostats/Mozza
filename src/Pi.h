#include <Rcpp.h>
using namespace Rcpp;

#ifndef PI_class
#define PI_class

class Pi_family {
  public:
    double a, b, c, alpha;
    Pi_family(double a_, double b_, double c_, double alpha_) :
      a(a_), b(b_), c(c_), alpha(alpha_) {};

    double operator()(double x1, double x2) {
      return pow( 1 + a*x1*x1 + 2*b*x1*x2 + c*x2*x2, alpha);
    }
};


#endif
