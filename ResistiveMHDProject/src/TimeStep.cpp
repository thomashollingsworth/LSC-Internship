#include <cmath>
#include "TimeStep.h"


void get_time_step(Grid& grid, const SimulationConfig& cfg){
    set_c_h(grid,cfg,grid.Prim);
    grid.dt= cfg.C *std::min(grid.dx,grid.dy)/(grid.c_h);
}