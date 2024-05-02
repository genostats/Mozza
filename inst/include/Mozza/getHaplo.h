#include <Rcpp.h>
#include "Mozza.h"

#ifndef _mozza_getHaplo_
#define _mozza_getHaplo_

namespace mozza {

template<typename T>
mosaic & getHaplo(T & zyg);

template<>
inline mosaic & getHaplo<mosaic>(mosaic & zyg) {
  return zyg;
};

template<>
inline mosaic & getHaplo<mosaic *>(mosaic * & zyg) {
  return *zyg;
};

} // end namespace

#endif
