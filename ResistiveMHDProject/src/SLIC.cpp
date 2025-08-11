#include "SLIC.h"
#include "StateVector.h"

void set_delta_x(Grid& grid, const SimulationConfig& cfg){
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    for(size_t i=1;i<nx+2*g-1;i++){
        for(size_t j=1;j<ny+2*g-1;j++){
            grid.Delta(i,j)=0.5*(grid.U(i+1,j)- grid.U(i-1,j));

    

}}}