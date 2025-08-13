#ifndef C_F_H
#define C_F_H


#include "Grid.h"
#include "SimulationConfig.h"

using PSV=PrimitiveStateVector;


void calc_cf_x(double gamma, double& output, const PrimitiveStateVector& p);
void set_cf_x(Grid& grid, const SimulationConfig& cfg,const Array2D<PSV>& prim_array);

void calc_cf_y(double gamma, double& output, const PSV& p);
void set_cf_y(Grid& grid, const SimulationConfig& cfg,const Array2D<PSV>& prim_array);

void calc_cf_z(double gamma, double& output, const PrimitiveStateVector& p);
void set_cf_z(Grid& grid, const SimulationConfig& cfg,const Array2D<PSV>& prim_array);

void set_c_h(Grid& grid, const SimulationConfig& cfg,const Array2D<PSV>& prim_array);
#endif