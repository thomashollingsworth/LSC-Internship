#include <cmath>
#include "TimeStep.h"
#include <iostream>


void get_time_step(Grid& grid, const SimulationConfig& cfg){
    set_c_h(grid,cfg,grid.Prim);
    grid.dt= cfg.C * std::fmin(grid.dx,grid.dy)/(grid.c_h);

}