#include <Rcpp.h>
#include "Mozza.h"

#ifndef _mozza_getZygote_
#define _mozza_getZygote_

namespace mozza {

template<typename T>
zygote & getZygote(T & zyg);

template<>
inline zygote & getZygote<zygote>(zygote & zyg) {
  return zyg;
};

template<>
inline zygote & getZygote<zygote *>(zygote * & zyg) {
  return *zyg;
};

/*
template<>
inline zygote & getZygote<Rcpp::internal::generic_proxy<19>>(Rcpp::internal::generic_proxy<19> & zyg) {
  Rcpp::XPtr<mozza::zygote> zyg1 = zyg;
  return *zyg1;
}

template<>
inline zygote & getZygote<Rcpp::internal::generic_proxy<19>>(Rcpp::internal::generic_proxy<19> & zyg) {
  return *((Rcpp::XPtr<zygote>) zyg);
};
*/

/*
inline zygote & getZygote(zygote & zig) {
  return zig;
}
*/

} // end namespace

#endif
