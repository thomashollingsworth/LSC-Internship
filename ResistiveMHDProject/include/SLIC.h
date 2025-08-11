#ifndef SLIC_H
#define SLIC_H

#include "Grid.h"
#include "SimulationConfig.h"

//Each step has its own dedicated function and then a single function collects all 
void set_delta(Grid& grid, const SimulationConfig& cfg);

void set_r(Grid& grid, const SimulationConfig& cfg);

void do_VanLeerLimiting(Grid& grid, const SimulationConfig& cfg);

void set_ubar(Grid& grid, const SimulationConfig& cfg);

void set_ubar_flux(Grid& grid, const SimulationConfig& cfg);

void set_ubarplus(Grid& grid, const SimulationConfig& cfg);

void do_SLIC_update(Grid& grid, const SimulationConfig& cfg);







#endif