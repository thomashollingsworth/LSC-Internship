#ifndef SLIC_H
#define SLIC_H

#include "Grid.h"
#include "SimulationConfig.h"

//Each step has its own dedicated function and then a single function collects all 
void set_delta_x(Grid& grid, const SimulationConfig& cfg);
void set_delta_y(Grid& grid, const SimulationConfig& cfg);

void set_r_x(Grid& grid, const SimulationConfig& cfg);
void set_r_y(Grid& grid, const SimulationConfig& cfg);

void do_VanLeerLimiting(Grid& grid, const SimulationConfig& cfg);

void set_ubar(Grid& grid, const SimulationConfig& cfg);

void set_ubar_xflux(Grid& grid, const SimulationConfig& cfg);

void set_ubar_yflux(Grid& grid, const SimulationConfig& cfg);

void set_ubar_plusx(Grid& grid, const SimulationConfig& cfg);

void set_ubar_plusy(Grid& grid, const SimulationConfig& cfg);

void do_SLIC_xupdate(Grid& grid, const SimulationConfig& cfg);

void do_SLIC_yupdate(Grid& grid, const SimulationConfig& cfg);







#endif