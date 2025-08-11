#ifndef C_F_H
#define C_F_H


#include "Grid.h"
#include "SimulationConfig.h"

void calc_cf_x(double gamma, const PrimitiveStateVector prim_cell);
void set_cf_x(Grid& grid, const SimulationConfig& cfg,const PSV& prim_array);

void calc_cf_y(double gamma, const PrimitiveStateVector prim_cell);
void set_cf_y(Grid& grid, const SimulationConfig& cfg,const PSV& prim_array);

void calc_cf_z(double gamma, const PrimitiveStateVector prim_cell);
void set_cf_z(Grid& grid, const SimulationConfig& cfg,const PSV& prim_array);

void set_c_h(Grid& grid, const SimulationConfig& cfg,const PSV& prim_array);
#endif