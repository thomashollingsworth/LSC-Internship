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
#include "BoundaryConditions.h"



int main() {
//Grid will hold all 2D arrays and data for intermediate values
Grid grid(100,100);
//Simulation config is used to control parameters of the test
SimulationConfig cfg;
cfg.save_interval=1;
cfg.plot_interval=1;
cfg.save_directory="TRIAL";

cfg.initialcondition= SimulationConfig::InitialCondition::BrioWu_x;

//Initialise function sets up grid using the simulation config settings

initialise(grid,cfg);
std::cout<<"Initialised"<<std::endl;
double time=0;
double iteration=1;

do{
    
    get_time_step(grid,cfg);
    save_to_file(grid,cfg,iteration,time);
    
    time+=grid.dt;
    iteration+=1;
    
    

    //do_half_psi_update(grid,cfg);
    
    //do_SLIC_xupdate(grid,cfg);
    
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    for(size_t i=1;i<nx+2*g-1;i++){
        for(size_t j=1;j<ny+2*g-1;j++){

            grid.uBarL(i,j)=grid.U(i,j);
            grid.uBarR(i,j)=grid.U(i,j);

        }}


    do_HLLD_x_update(grid,cfg);

    grid.reset_intermediate_arrays();

    //do_SLIC_yupdate(grid,cfg);

    //do_HLLD_y_update(grid,cfg);

    //grid.reset_intermediate_arrays();

    //do_half_psi_update(grid,cfg);

    

    std::cout<<"Iteration " <<iteration << " completed t= "<<time<<std::endl;
}while(time<cfg.t_end);

}
