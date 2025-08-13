#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

#include "StateVector.h"
#include "Grid.h"
#include "SimulationConfig.h"

template<typename T>
void update_bcs(Grid& grid, const SimulationConfig& cfg, Array2D<T>& array){
    
    
    using BCs=SimulationConfig::BoundaryCondition;
    const size_t nx = grid.num_xcells;
    const size_t ny = grid.num_ycells;
    const size_t g  = grid.ghost_cells;
    //Update each edge in turn
   
    switch(cfg.bcs_x0){
        case BCs::Periodic:
        for(size_t j=0;j<(ny+2*g);j++){
            array(0,j)=array(nx,j);
            array(1,j)=array(nx+1,j);
        }
        break;
        case BCs::Transmissive:
        for(size_t j=0;j<(ny+2*g);j++){
            array(0,j)=array(g,j);
            array(1,j)=array(g,j);
        }
        break;
        case BCs::Reflective:
        for(size_t j=0;j<(ny+2*g);j++){
            array(0,j)=array(3,j);
            array(1,j)=array(2,j);
            if constexpr (std::is_same_v<T, ConservedStateVector>) {
            //Reflect the normal vector quantities
            array(0,j).B().x()*=-1;
            array(1,j).B().x()*=-1;
            array(0,j).momentum().x()*=-1;
            array(1,j).momentum().x()*=-1;}
            else if constexpr (std::is_same_v<T, PrimitiveStateVector>){
            array(0,j).B().x()*=-1;
            array(1,j).B().x()*=-1;
            array(0,j).velocity().x()*=-1;
            array(1,j).velocity().x()*=-1;}
            }
            break;
        }

    switch (cfg.bcs_xf) {
        case BCs::Periodic:
            for (size_t j = 0; j < ny + 2*g; j++) {
                array(nx + g,     j) = array(g, j);
                array(nx + g + 1, j) = array(g + 1, j);
            }
            break;
        case BCs::Transmissive:
            for (size_t j = 0; j < ny + 2*g; j++) {
                array(nx + g,     j) = array(nx + g - 1, j);
                array(nx + g + 1, j) = array(nx + g - 1, j);
            }
            break;
        case BCs::Reflective:
            for (size_t j = 0; j < ny + 2*g; j++) {
                array(nx + g,     j) = array(nx + g - 1, j);
                array(nx + g + 1, j) = array(nx + g - 2, j);
                if constexpr (std::is_same_v<T, ConservedStateVector>) {
                //Reflect the normal vector quantities
                array(nx+g,j).B().x()*=-1;
                array(nx+g+1,j).B().x()*=-1;
                array(nx+g,j).momentum().x()*=-1;
                array(nx+g+1,j).momentum().x()*=-1;}
                else if constexpr (std::is_same_v<T, PrimitiveStateVector>){
                array(nx+g,j).B().x()*=-1;
                array(nx+g+1,j).B().x()*=-1;
                array(nx+g,j).velocity().x()*=-1;
                array(nx+g+1,j).velocity().x()*=-1;}
            }
            
            break;
    }
    
    
    switch(cfg.bcs_y0){
        case BCs::Periodic:
        for(size_t i=0;i<(nx+2*g);i++){
            array(i,0)=array(i,ny);
            array(i,1)=array(i,ny+1);
        }
        break;
        case BCs::Transmissive:
        for(size_t i=0;i<(nx+2*g);i++){
            array(i,0)=array(i,g);
            array(i,1)=array(i,g);
        }
        break;
        case BCs::Reflective:
        for(size_t i=0;i<(nx+2*g);i++){
            array(i,0)=array(i,g+1);
            array(i,1)=array(i,g);
            if constexpr (std::is_same_v<T, ConservedStateVector>) {
            //Reflect the normal vector quantities
            array(i,0).B().y()*=-1;
            array(i,1).B().y()*=-1;
            array(i,0).momentum().y()*=-1;
            array(i,1).momentum().y()*=-1;}
            else if constexpr (std::is_same_v<T, PrimitiveStateVector>){
            array(i,0).B().y()*=-1;
            array(i,1).B().y()*=-1;
            array(i,0).velocity().y()*=-1;
            array(i,1).velocity().y()*=-1;}
            }
            break;
        }

    switch (cfg.bcs_yf) {
        case BCs::Periodic:
            for (size_t i = 0; i < nx + 2*g; i++) {
                array(i,ny+g) = array(i,g);
                array(i,ny + g + 1) = array(i,g+1);
            }
            break;
        case BCs::Transmissive:
            for (size_t i = 0; i < nx + 2*g; i++) {
                array(i,ny+g) = array(i,ny+g-1);
                array(i,ny+g+1) = array(i,ny+g-1);
            }
            break;
        case BCs::Reflective:
            for (size_t i = 0; i < nx + 2*g; i++) {
                array(i,ny+g) = array(i,ny+g-1);
                array(i,ny+g+1) = array(i,ny+g-2);
                if constexpr (std::is_same_v<T, ConservedStateVector>) {
                //Reflect the normal vector quantities
                array(i,ny+g).B().y()*=-1;
                array(i,ny+g+1).B().y()*=-1;
                array(i,ny+g).momentum().y()*=-1;
                array(i,ny+g+1).momentum().y()*=-1;}
                else if constexpr (std::is_same_v<T, PrimitiveStateVector>){
                array(i,ny+g).B().y()*=-1;
                array(i,ny+g+1).B().y()*=-1;
                array(i,ny+g).velocity().y()*=-1;
                array(i,ny+g+1).velocity().y()*=-1;}
            }
            break;
    }
}

#endif