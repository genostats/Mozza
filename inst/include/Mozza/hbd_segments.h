#include <Rcpp.h>
#include "mosaic.h"
#include "zygote.h"
#include "HBD_at_point.h"
#include "segments.h"
/* 
 * calcul des segments HBD
 */

#ifndef _mozza_hbd_segments_
#define _mozza_hbd_segments_

namespace mozza {

// segments partagés HBD par les deux haplotypes du zygote
segments HBD_segments(zygote & Z);

}

#endif
