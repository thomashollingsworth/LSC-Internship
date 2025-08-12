#include "psi_source_term.h"
#include <cmath>

void do_half_psi_update(Grid& grid,const SimulationConfig& cfg){
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    double c_p_squared=grid.c_h*0.18;
    
    for(size_t i=0;i<nx+2*g;i++){
        for(size_t j=0;j<ny+2*g;j++){
            grid.U(i,j).psi()*=std::exp((-grid.dt/2)*(grid.c_h*grid.c_h)/c_p_squared);
            grid.Prim(i,j).psi()*=std::exp((-grid.dt/2)*(grid.c_h*grid.c_h)/c_p_squared);

}}
}
