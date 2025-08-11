#ifndef SAVINGROUTINE_H
#define SAVINGROUTINE_H

#include <filesystem>  // C++17
#include <fstream>
#include <iostream>

#include "Grid.h"
#include "SimulationConfig.h"

//Run at the end of every step, checks whether  
void save_to_file(const Grid& grid, const SimulationConfig& cfg, size_t iteration,double time);

#endif