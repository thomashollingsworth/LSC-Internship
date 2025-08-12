#include "Grid.h"
#include "SimulationConfig.h"
#include "Initialise.h"
#include "SavingRoutine.h"
#include "MHD_flux.h"
#include "c_f.h"
#include "TimeStep.h"
#include "SLIC.h"
#include "HLLD.h"
#include "psi_source_term.h"
#include <iostream>



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

do{
    get_time_step(grid,cfg);
    time+=grid.dt;
    iteration+=1;
    save_to_file(grid,cfg,iteration,time);

    do_half_psi_update(grid,cfg);
    
    do_SLIC_xupdate(grid,cfg);
    do_HLLD_x_update(grid,cfg);
    
    grid.reset_intermediate_arrays();

    do_SLIC_yupdate(grid,cfg);
    do_HLLD_y_update(grid,cfg);

    do_half_psi_update(grid,cfg);

    std::cout<<"Iteration " <<iteration << "completed"<<std::endl;
}while(time<cfg.t_end);

}
