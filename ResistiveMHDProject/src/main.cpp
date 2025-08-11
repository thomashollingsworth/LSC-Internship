#include "Grid.h"
#include "SimulationConfig.h"
#include "Initialise.h"
#include "SavingRoutine.h"
#include "MHD_flux.h"


int main() {
//Grid will hold all 2D arrays and data for intermediate values
Grid grid(256,256);
//Simulation config is used to control parameters of the test
SimulationConfig cfg;
cfg.save=true;
cfg.save_directory="TRIAL";
cfg.plot=true;

cfg.initialcondition= SimulationConfig::InitialCondition::OrzsagTang;

//Initialise function sets up grid using the simulation config settings

initialise(grid,cfg);
std::cout<<"initialised"<<std::endl;
double time=0;
double iteration=0;
save_to_file(grid,cfg,iteration,time);

}
