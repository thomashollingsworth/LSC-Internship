#ifndef SIMULATIONCONFIG_H
#define SIMULATIONCONFIG_H

#include <string>
//Structure that is used to choose the initial conditions, boundary conditions etc for the simulation

struct SimulationConfig{

    SimulationConfig(){};
    
    enum class InitialCondition{
        BrioWu_x,
        BrioWu_y,
        OrzsagTang,
        KelvinHelmholtz,
    };

    enum class BoundaryCondition{
        Transimssive,
        Periodic,
        Reflective
    };

    //Test parameter members initialised with default values

    InitialCondition initialcondition=InitialCondition::OrzsagTang;
    double t_end=1.0;
    double C=0.8;
    double gamma=5./3.;
    BoundaryCondition bcs_x0=BoundaryCondition::Periodic;
    BoundaryCondition bcs_xf=BoundaryCondition::Periodic;
    BoundaryCondition bcs_y0=BoundaryCondition::Periodic;
    BoundaryCondition bcs_yf=BoundaryCondition::Periodic;

    //Saving Parameters set to a default
    bool save=true;
    int save_interval=10;
    bool plot=true;
    int plot_interval=20;

    std::string save_directory="Results";

};











#endif