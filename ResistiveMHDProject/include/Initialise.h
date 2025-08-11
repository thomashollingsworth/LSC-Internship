#ifndef INITIALISE_H
#define INITIALISE_H


#include "Grid.h"
#include "SimulationConfig.h"
#include "BoundaryConditions.h"


//Function to be called once config settings and grid dimensions have been decided

void initialise(Grid& grid, SimulationConfig& cfg);

#endif