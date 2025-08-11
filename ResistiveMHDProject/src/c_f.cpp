#include "c_f.h"
#include <cmath>

using PSV=Array2D<PrimitiveStateVector>;

void calc_cf_x(double gamma, double& output, const PrimitiveStateVector& p){
    double factor= (gamma*p.pressure()+dot(p.B(),p.B()))/p.density();
    output=std::sqrt(0.5*(factor+std::sqrt(factor*factor-4.0*gamma*p.pressure()*p.B().x()*p.B().x()/(p.density()*p.density()))));
}

void set_cf_x(Grid& grid, const SimulationConfig& cfg,const PSV& prim_array){
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    for(size_t i=0;i<nx+2*g;i++){
        for(size_t j=0;j<ny+2*g;j++){
            calc_cf_x(cfg.gamma,grid.c_fx(i,j),prim_array(i,j));
    }};
}


void calc_cf_y(double gamma, double& output, const PrimitiveStateVector& p){
    double factor= (gamma*p.pressure()+dot(p.B(),p.B()))/p.density();
    output=std::sqrt(0.5*(factor+std::sqrt(factor*factor-4.0*gamma*p.pressure()*p.B().y()*p.B().y()/(p.density()*p.density()))));
}

void set_cf_y(Grid& grid, const SimulationConfig& cfg,const PSV& prim_array){
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    for(size_t i=0;i<nx+2*g;i++){
        for(size_t j=0;j<ny+2*g;j++){
            calc_cf_y(cfg.gamma,grid.c_fy(i,j),prim_array(i,j));
    }};
}


void calc_cf_z(double gamma, double& output, const PrimitiveStateVector& p){
    double factor= (gamma*p.pressure()+dot(p.B(),p.B()))/p.density();
    output=std::sqrt(0.5*(factor+std::sqrt(factor*factor-4.0*gamma*p.pressure()*p.B().z()*p.B().z()/(p.density()*p.density()))));
}

void set_cf_z(Grid& grid, const SimulationConfig& cfg,const PSV& prim_array){

    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    for(size_t i=0;i<nx+2*g;i++){
        for(size_t j=0;j<ny+2*g;j++){
            calc_cf_z(cfg.gamma,grid.c_fz(i,j),prim_array(i,j));
    }};
}

void set_c_h(Grid& grid, const SimulationConfig& cfg,const PSV& prim_array){
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    set_cf_x(grid,cfg,prim_array);
    set_cf_y(grid,cfg,prim_array);
    set_cf_z(grid,cfg,prim_array);
    
    double max_ch=0.;
    
    for(size_t i=0;i<nx+2*g;i++){
        for(size_t j=0;j<ny+2*g;j++){

        double new_max=std::max({
            std::abs(prim_array(i,j).velocity().x())+grid.c_fx(i,j),
            std::abs(prim_array(i,j).velocity().y())+grid.c_fy(i,j),
            std::abs(prim_array(i,j).velocity().z())+grid.c_fz(i,j)});

            if(new_max> max_ch){
            max_ch=new_max;
        } 



        }}

        grid.c_h= max_ch;
}

