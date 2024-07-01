#include "mosaic.h"

#ifndef ZYGOTE__
#define ZYGOTE__
namespace mozza {
class zygote : public std::pair<mosaic, mosaic> {
public:
  // constructeur avec deux haplotypes mosaic
  zygote(mosaic h1, mosaic h2) : std::pair<mosaic, mosaic>(h1, h2)  {}
  // deux constructeur avec les paramètres de constructeurs de mosaic
  zygote(const std::vector<double> & chr_len, int ntiles, double mean_length)
    : std::pair<mosaic, mosaic>(mosaic(chr_len, ntiles, mean_length), mosaic(chr_len, ntiles, mean_length)) {}
  template<typename vec>
  zygote(const std::vector<double> & chr_len, const vec & proba_tiles, double mean_length)
    : std::pair<mosaic, mosaic>(mosaic(chr_len, proba_tiles, mean_length), mosaic(chr_len, proba_tiles, mean_length)) {}
  // un constructeur avec deux zygotes : reproduction
  zygote(zygote & Z1, zygote & Z2) : std::pair<mosaic, mosaic>(mosaic(Z1.first, Z1.second), mosaic(Z2.first, Z2.second)) {}
};

inline zygote operator+(zygote Z1, zygote Z2) {
  return zygote(mosaic(Z1.first, Z1.second), mosaic(Z2.first, Z2.second));
}

}
#endif
