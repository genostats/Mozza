#include "mosaic.h"

#ifndef ZYGOTE__
#define ZYGOTE__
namespace mozza {
class zygote : public std::pair<mosaic, mosaic> {
public:
  zygote(mosaic h1, mosaic h2) : std::pair<mosaic, mosaic>(h1, h2)  {}
  zygote(const std::vector<double> & chr_len, int ntiles, double mean_length) :
  std::pair<mosaic, mosaic>(mosaic(chr_len, ntiles, mean_length), mosaic(chr_len, ntiles, mean_length)) {}
};

inline zygote operator+(zygote Z1, zygote Z2) {
  return zygote(mosaic(Z1.first, Z1.second), mosaic(Z2.first, Z2.second));
}

};
#endif
