#include "SLIC.h"
#include "StateVector.h"
#include "MHD_flux.h"
#include "BoundaryConditions.h"

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
    //Slope limiting using primitive variables 
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    for(size_t i=1;i<nx+2*g-1;i++){
        for(size_t j=1;j<ny+2*g-1;j++){ 
            
            grid.Chi(i,j)=StateVector((grid.Prim(i,j)- grid.Prim(i-1,j))/((grid.Prim(i+1,j)- grid.Prim(i,j))+1e-12));
           
}}}

void set_r_y(Grid& grid, const SimulationConfig& cfg){
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    for(size_t i=1;i<nx+2*g-1;i++){
        for(size_t j=1;j<ny+2*g-1;j++){
            grid.Chi(i,j)=StateVector((grid.Prim(i,j)- grid.Prim(i,j-1))/((grid.Prim(i,j+1)- grid.Prim(i,j))+1e-12));
            
            
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
            grid.uBarL(i,j)+=grid.U(i,j);
            grid.uBarR(i,j)+=grid.U(i,j);
            grid.uBarL(i,j)-=0.5*grid.Chi(i,j)*grid.Delta(i,j);
            grid.uBarR(i,j)+=0.5*grid.Chi(i,j)*grid.Delta(i,j);
}}
}

void set_ubar_xflux(Grid& grid, const SimulationConfig& cfg){
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    for(size_t i=1;i<nx+2*g-1;i++){
        for(size_t j=1;j<ny+2*g-1;j++){
        //Set the flux L&R arrays
        MHD_xflux(grid.fluxL(i,j),grid.uBarL(i,j),grid.uBarL(i,j).con_to_prim(cfg.gamma),grid.c_h);
        MHD_xflux(grid.fluxR(i,j),grid.uBarR(i,j),grid.uBarR(i,j).con_to_prim(cfg.gamma),grid.c_h);
    }}
}

void set_ubar_yflux(Grid& grid, const SimulationConfig& cfg){
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    for(size_t i=1;i<nx+2*g-1;i++){
        for(size_t j=1;j<ny+2*g-1;j++){
        //Set the flux L&R arrays
        MHD_yflux(grid.fluxL(i,j),grid.uBarL(i,j),grid.uBarL(i,j).con_to_prim(cfg.gamma),grid.c_h);
        MHD_yflux(grid.fluxR(i,j),grid.uBarR(i,j),grid.uBarR(i,j).con_to_prim(cfg.gamma),grid.c_h);
    }}
}

void set_ubar_plusx(Grid& grid, const SimulationConfig& cfg){
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    for(size_t i=1;i<nx+2*g-1;i++){
        for(size_t j=1;j<ny+2*g-1;j++){
            grid.uBarL(i,j)-=0.5*grid.dt/grid.dx*(grid.fluxR(i,j)-grid.fluxL(i,j));
            grid.uBarR(i,j)-=0.5*grid.dt/grid.dx*(grid.fluxR(i,j)-grid.fluxL(i,j));
}}
update_bcs(grid,cfg,grid.uBarL);
update_bcs(grid,cfg,grid.uBarR);
}

void set_ubar_plusy(Grid& grid, const SimulationConfig& cfg){
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    for(size_t i=1;i<nx+2*g-1;i++){
        for(size_t j=1;j<ny+2*g-1;j++){
            grid.uBarL(i,j)-=0.5*grid.dt/grid.dy*(grid.fluxR(i,j)-grid.fluxL(i,j));
            grid.uBarR(i,j)-=0.5*grid.dt/grid.dy*(grid.fluxR(i,j)-grid.fluxL(i,j));
}}
update_bcs(grid,cfg,grid.uBarL);
update_bcs(grid,cfg,grid.uBarR);
}

void do_SLIC_xupdate(Grid& grid, const SimulationConfig& cfg){
    //Performs entire SLIC process starting from U and ending with 
    //putting Ubarplus arrays in place
    set_delta_x(grid,cfg);
    set_r_x(grid,cfg);
    do_VanLeerLimiting(grid,cfg);
    set_ubar(grid,cfg);
    set_ubar_xflux(grid,cfg);
    set_ubar_plusx(grid,cfg);
}
void do_SLIC_xupdate(Grid& grid, const SimulationConfig& cfg){
    //Performs entire SLIC process starting from U and ending with -
    //putting Ubarplus arrays in place
    set_delta_y(grid,cfg);
    set_r_y(grid,cfg);
    do_VanLeerLimiting(grid,cfg);
    set_ubar(grid,cfg);
    set_ubar_yflux(grid,cfg);
    set_ubar_plusy(grid,cfg);
}