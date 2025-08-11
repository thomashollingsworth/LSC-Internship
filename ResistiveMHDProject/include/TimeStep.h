#ifndef TIMESTEP_H
#define TIMESTEP_H
#include "c_f.h"
#include "Grid.h"
#include "SimulationConfig.h"

void get_time_step(Grid& grid, const SimulationConfig& cfg);


#endif