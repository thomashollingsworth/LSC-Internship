#include "SLIC.h"
#include "StateVector.h"

void set_delta_x(Grid& grid, const SimulationConfig& cfg){
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    for(size_t i=1;i<nx+2*g-1;i++){
        for(size_t j=1;j<ny+2*g-1;j++){
            grid.Delta(i,j)=StateVector(0.5*(grid.U(i+1,j)- grid.U(i-1,j)));
}}}

void set_delta_y(Grid& grid, const SimulationConfig& cfg){
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    for(size_t i=1;i<nx+2*g-1;i++){
        for(size_t j=1;j<ny+2*g-1;j++){
            grid.Delta(i,j)=StateVector(0.5*(grid.U(i,j+1)- grid.U(i,j-1)));
}}}


void set_r_x(Grid& grid, const SimulationConfig& cfg){
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    for(size_t i=1;i<nx+2*g-1;i++){
        for(size_t j=1;j<ny+2*g-1;j++){ 
            
            grid.Chi(i,j)=StateVector((grid.U(i,j)- grid.U(i-1,j))/((grid.U(i+1,j)- grid.U(i,j))+1e-12));
            
            
}}}

void set_r_y(Grid& grid, const SimulationConfig& cfg){
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    for(size_t i=1;i<nx+2*g-1;i++){
        for(size_t j=1;j<ny+2*g-1;j++){
            grid.Chi(i,j)=StateVector((grid.U(i,j)- grid.U(i,j-1))/((grid.U(i,j+1)- grid.U(i,j))+1e-12));
            
            
}}}

void do_VanLeerLimiting(Grid& grid, const SimulationConfig& cfg){
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    for(size_t i=1;i<nx+2*g-1;i++){
        for(size_t j=1;j<ny+2*g-1;j++){
            for(size_t k=0;j<9;j++){
                if(grid.Chi(i,j)[k]<=0){
                grid.Chi(i,j)[k]*=0;} 
                else{
                grid.Chi(i,j)[k]=std::min(2*grid.Chi(i,j)[k]/(1+grid.Chi(i,j)[k]),2/(1+grid.Chi(i,j)[k]));} 
    }}}
}

void set_ubar(Grid& grid, const SimulationConfig& cfg){
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    for(size_t i=1;i<nx+2*g-1;i++){
        for(size_t j=1;j<ny+2*g-1;j++){
            grid.uBarL(i,j)*=0;//reset ubar arrays
            grid.uBarR(i,j)*=0;
            grid.uBarL(i,j)+=grid.U(i,j);
            grid.uBarL(i,j)-=0.5*grid.Chi(i,j)*grid.Delta(i,j);
            grid.uBarR(i,j)+=0.5*grid.Chi(i,j)*grid.Delta(i,j);





}}
