#ifndef GRID_H
#define GRID_H


#include "StateVector.h"
#include "Vector3.h"
#include "Array2D.h"
#include <vector>


//Grid structure: 
//Takes num_x and num_y as constructor parameters
//Has all persistent and intermediate arrays as members
//Initialises arrays on construction 
//Contains a template function for converting between 2D indices and 1D arrays
//Contains various template functions for boundary conditions

struct Grid{
    
    size_t ghost_cells; //Won't realistically be changing 
    size_t num_xcells,num_ycells;

    //Constructor

    Grid(size_t nx, size_t ny,size_t g_=2);

    //These values will be calculated in initialise function
    std::vector<double> x;
    std::vector<double> y;

    double dx;
    double dy;
    
    //Persistent arrays
    Array2D<ConservedStateVector> U;
    Array2D<PrimitiveStateVector> Prim;
    
    

    //Intermediate arrays (memory buffers)

    //  - FOR SLOPE LIMITING
    Array2D<ConservedStateVector> uBarL;
    Array2D<ConservedStateVector> uBarR;

    Array2D<StateVector> Delta;
    Array2D<StateVector> Chi;
    
    Array2D<StateVector> fluxL;
    Array2D<StateVector> fluxR;

     //  - FOR MUSCL-HANCOCK SCHEME
    Array2D<double> c_fx;
    Array2D<double> c_fy;
    Array2D<double> c_fz;
    
    Array2D<double> S_L;
    Array2D<double> S_R;
    Array2D<double> S_M;

    double c_h;
    double dt;
};


















#endif