#ifndef HLLD_H
#define HLLD_H

#include "Grid.h"
#include "SimulationConfig.h"
#include "c_f.h"

//used for deciding what flux each point should take
enum class Region{
    L,
    L_star,
    L_doublestar,
    R_doublestar,
    R_star,
    R
};

//Update the uBar arrays in placce
void set_xtilde_vals(Grid& grid, const SimulationConfig& cfg);
void set_ytilde_vals(Grid& grid, const SimulationConfig& cfg);

void update_fluxes_x(Grid& grid, const SimulationConfig& cfg);
void update_fluxes_y(Grid& grid, const SimulationConfig& cfg);

//Update the primL/R arrays to match the uBar arrays
void set_prim_bar(Grid& grid, const SimulationConfig& cfg);




//Single cell level operations (storing doubles is cheap and easy)

std::tuple<double,double> calc_S_LR_x(const PrimitiveStateVector& pL,const PrimitiveStateVector& pR,const double gamma);
//(use calc_cfx for primL and primR)
std::tuple<double,double> calc_S_LR_y(const PrimitiveStateVector& pL,const PrimitiveStateVector& pR,const double gamma);
//(use calc_cfx for primL and primR)


double calc_S_M_x(const double S_L,const double S_R,const PrimitiveStateVector& pL,const PrimitiveStateVector& pR);
double calc_S_M_y(const double S_L,const double S_R,const PrimitiveStateVector& pL,const PrimitiveStateVector& pR);

std::tuple<double,double> calc_S_stars_x(const double S_L,const double S_R,const double S_M,const PrimitiveStateVector& pL,const PrimitiveStateVector& pR,const double gamma);
//Requires calc of density star states 
std::tuple<double,double> calc_S_stars_y(const double S_L,const double S_R,const double S_M,const PrimitiveStateVector& pL,const PrimitiveStateVector& pR,const double gamma);


void calc_u_starLR_x(bool L,ConservedStateVector& output,const double S_M,const double S_L,const double S_R,const PrimitiveStateVector& pL,const PrimitiveStateVector& pR,double gamma);
void calc_u_starLR_y(bool L,ConservedStateVector& output,const double S_M,const double S_L,const double S_R,const PrimitiveStateVector& pL,const PrimitiveStateVector& pR,double gamma);
//Full in place star state calc, can be either L or R state
//Must do L or R for star regions, L and R for double star regions

void calc_u_doublestarLR_x(bool L,ConservedStateVector& u_starL,ConservedStateVector& u_starR);
void calc_u_doublestarLR_y(bool L,ConservedStateVector& u_starL,ConservedStateVector& u_starR);

void calc_HLLD_flux_x(StateVector& flux_out,const StateVector& flux_L,const StateVector& flux_R, ConservedStateVector& u_starL,ConservedStateVector& u_starR,ConservedStateVector& u_doublestar,const PrimitiveStateVector& pL,const PrimitiveStateVector& pR,const ConservedStateVector& uL,const ConservedStateVector& uR,const double gamma);

void do_HLLD_x_update(Grid& grid, const SimulationConfig& cfg);
void do_HLLD_y_update(Grid& grid, const SimulationConfig& cfg);





#endif