#ifndef MHD_FLUX_H
#define MHD_FLUX_H

#include "StateVector.h"

void MHD_xflux(StateVector& out, const ConservedStateVector& CSV, const PrimitiveStateVector& PSV,const double c_h);
void MHD_yflux(StateVector& out, const ConservedStateVector& CSV, const PrimitiveStateVector& PSV,const double c_h);

#endif