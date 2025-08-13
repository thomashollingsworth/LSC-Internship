#include "HLLD.h"
#include <cmath>
#include "MHD_flux.h"
#include "BoundaryConditions.h"

//L and R states in HLLD are defined differently to those in SLIC
void redefine_LRx(Grid& grid,const SimulationConfig& cfg){
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    for(size_t i=1;i<nx+2*g-2;i++){
        for(size_t j=1;j<ny+2*g-2;j++){
            grid.primL(i,j)=grid.uBarR(i,j).con_to_prim(cfg.gamma);
            grid.primR(i,j)=grid.uBarL(i+1,j).con_to_prim(cfg.gamma);

}}
}

void redefine_LRy(Grid& grid,const SimulationConfig& cfg){
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    for(size_t i=1;i<nx+2*g-2;i++){
        for(size_t j=1;j<ny+2*g-2;j++){
            grid.primL(i,j)=grid.uBarR(i,j).con_to_prim(cfg.gamma);
            grid.primR(i,j)=grid.uBarL(i,j+1).con_to_prim(cfg.gamma);

}}
}

void set_xtilde_vals(Grid& grid){
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    for(size_t i=1;i<nx+2*g-2;i++){
        for(size_t j=1;j<ny+2*g-2;j++){
            grid.primL(i,j).B().x()=0.5*(grid.primL(i,j).B().x()+grid.primR(i,j).B().x())-0.5*1./grid.c_h*(grid.primR(i,j).psi()-grid.primL(i,j).psi());
            grid.primL(i,j).psi()=0.5*(grid.primR(i,j).psi()+grid.primL(i,j).psi())-0.5*grid.c_h*(grid.primR(i,j).B().x()-grid.primL(i,j).B().x());
            grid.primR(i,j).B().x()=grid.primL(i,j).B().x();
            grid.primR(i,j).psi()=grid.primL(i,j).psi();

}}
}
void set_ytilde_vals(Grid& grid){
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    for(size_t i=1;i<nx+2*g-2;i++){
        for(size_t j=1;j<ny+2*g-2;j++){
            grid.primL(i,j).B().y()=0.5*(grid.primL(i,j).B().y()+grid.primR(i,j).B().y())-0.5*1./grid.c_h*(grid.primR(i,j).psi()-grid.primL(i,j).psi());
            grid.primL(i,j).psi()=0.5*(grid.primR(i,j).psi()+grid.primL(i,j).psi())-0.5*grid.c_h*(grid.primR(i,j).B().y()-grid.primL(i,j).B().y());
            grid.primR(i,j).B().y()=grid.primL(i,j).B().y();
            grid.primR(i,j).psi()=grid.primL(i,j).psi();

}}
}

void set_u_bar(Grid& grid, const SimulationConfig& cfg){
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    
    for(size_t i=1;i<nx+2*g-2;i++){
        for(size_t j=1;j<ny+2*g-2;j++){
            grid.uBarL(i,j)=grid.primL(i,j).prim_to_con(cfg.gamma);
            grid.uBarR(i,j)=grid.primR(i,j).prim_to_con(cfg.gamma);
}}
}

void update_fluxes_x(Grid& grid){
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    for(size_t i=1;i<nx+2*g-2;i++){
        for(size_t j=1;j<ny+2*g-2;j++){
            MHD_xflux(grid.fluxR(i,j),grid.uBarR(i,j),grid.primR(i,j),grid.c_h);
            MHD_xflux(grid.fluxL(i,j),grid.uBarL(i,j),grid.primL(i,j),grid.c_h);
}}
}

void update_fluxes_y(Grid& grid){
    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    for(size_t i=1;i<nx+2*g-2;i++){
        for(size_t j=1;j<ny+2*g-2;j++){
            MHD_yflux(grid.fluxR(i,j),grid.uBarR(i,j),grid.primR(i,j),grid.c_h);
            MHD_yflux(grid.fluxL(i,j),grid.uBarL(i,j),grid.primL(i,j),grid.c_h);
}}
}

std::tuple<double,double> calc_S_LR_x(const PrimitiveStateVector& pL,const PrimitiveStateVector& pR,const double gamma){
    double c_fL;
    double c_fR;
    calc_cf_x(gamma, c_fL, pL);
    calc_cf_x(gamma, c_fR, pR);
    double SL=std::min(pL.velocity().x(),pR.velocity().x())-std::max(c_fL,c_fR);
   double SR=std::max(pL.velocity().x(),pR.velocity().x())-std::min(c_fL,c_fR);
   return{SL,SR};
}

std::tuple<double,double> calc_S_LR_y(const PrimitiveStateVector& pL,const PrimitiveStateVector& pR,const double gamma){
    double c_fL;
    double c_fR;
    calc_cf_y(gamma, c_fL, pL);
    calc_cf_y(gamma, c_fR, pR);
    double SL=std::min(pL.velocity().y(),pR.velocity().y())-std::max(c_fL,c_fR);
   double SR=std::max(pL.velocity().y(),pR.velocity().y())-std::min(c_fL,c_fR);
   return{SL,SR};
}


double calc_S_M_x(const double S_L,const double S_R,const PrimitiveStateVector& pL,const PrimitiveStateVector& pR){
    double S_M;
    S_M= ((S_R-pR.velocity().x())*pR.density()*pR.velocity().x()-(S_L-pL.velocity().x())*pL.density()*pL.velocity().x() -pR.pressure_T()+pL.pressure_T())/((S_R-pR.velocity().x())*pR.density()-(S_L-pL.velocity().x())*pL.density());
    return S_M;
}

double calc_S_M_y(const double S_L,const double S_R,const PrimitiveStateVector& pL,const PrimitiveStateVector& pR){
    double S_M;
    S_M= ((S_R-pR.velocity().y())*pR.density()*pR.velocity().y()-(S_L-pL.velocity().x())*pL.density()*pL.velocity().y() -pR.pressure_T()+pL.pressure_T())/((S_R-pR.velocity().y())*pR.density()-(S_L-pL.velocity().y())*pL.density());
    return S_M;
}

std::tuple<double,double> calc_S_stars_x(const double S_L,const double S_R,const double S_M,const PrimitiveStateVector& pL,const PrimitiveStateVector& pR){
//Requires calc of density star states
double density_starL=pL.density()*(S_L-pL.velocity().x())/(S_L-S_M);
double density_starR=pR.density()*(S_R-pR.velocity().x())/(S_R-S_M);

double S_starL= S_M-std::fabs(pL.B().x())/std::sqrt(density_starL);
double S_starR= S_M+std::fabs(pR.B().x())/std::sqrt(density_starR);

return {S_starL,S_starR};
}

std::tuple<double,double> calc_S_stars_y(const double S_L,const double S_R,const double S_M,const PrimitiveStateVector& pL,const PrimitiveStateVector& pR){
//Requires calc of density star states
double density_starL=pL.density()*(S_L-pL.velocity().y())/(S_L-S_M);
double density_starR=pR.density()*(S_R-pR.velocity().y())/(S_R-S_M);

double S_starL= S_M-std::fabs(pL.B().y())/std::sqrt(density_starL);
double S_starR= S_M+std::fabs(pR.B().y())/std::sqrt(density_starR);

return {S_starL,S_starR};
}

void calc_u_starLR_x(bool L,ConservedStateVector& output,const double S_M,const double S_L,const double S_R,const PrimitiveStateVector& pL,const PrimitiveStateVector& pR,double gamma){
    double TOL=1e-8; //for detecting 0 denominator
    const PrimitiveStateVector& p = (L ? pL : pR);
    double S = (L ? S_L : S_R);

    output.density()=p.density()*(S-p.velocity().x())/(S-S_M);
    output.momentum().x()=output.density()*S_M;
    //Check denominators
    double vel_denom= (p.density()*(S-p.velocity().x())*(S-S_M)-p.B().x()*p.B().x());
    if (vel_denom<TOL){
        output.momentum().y()=output.density()*p.velocity().y();
        output.momentum().z()=output.density()*p.velocity().z();
    }
    else{
    output.momentum().y()= output.density()*(p.velocity().y()-((S_M-p.velocity().x())*p.B().x()*p.B().y())/vel_denom);
    output.momentum().z()= output.density()*(p.velocity().z()-((S_M-p.velocity().x())*p.B().x()*p.B().z())/vel_denom);
    }
    output.B().x()=p.B().x();

    double B_denom=(p.density()*(S-p.velocity().x())*(S-S_M)-p.B().x()*p.B().x());
    if (B_denom<TOL){
        output.B().y()=0.;
        output.B().z()=0.;
    }
    else{
    output.B().y()=p.B().y()*(p.density()*(S-p.velocity().x())*(S-p.velocity().x())-p.B().x()*p.B().x())/B_denom;
    output.B().z()=p.B().z()*(p.density()*(S-p.velocity().x())*(S-p.velocity().x())-p.B().x()*p.B().x())/B_denom;
    }
    
    double p_T_star=((S_R-pR.velocity().x())*pR.density()*pL.pressure_T()-(S_L-pL.velocity().x())*pL.density()*pR.pressure_T()+pL.density()*pR.density()*(S_R-pR.velocity().x())*(S_L-pL.velocity().x())*(pR.velocity().x()-pL.velocity().x()))/((S_R-pR.velocity().x())*pR.density()-(S_L-pL.velocity().x())*pL.density());
    //Intermediate step for convenience
    double E=p.pressure()/(gamma-1)+0.5*p.density()*dot(p.velocity(),p.velocity())+ 0.5* dot(p.B(),p.B());

    output.energy()=((S-p.velocity().x())*E-p.pressure_T()*p.velocity().x()+p_T_star*S_M + p.B().x()*(dot(p.velocity(),p.B())-1/output.density()*dot(output.momentum(),output.B())))/(S-S_M);
    output.psi()=p.psi();

}

void calc_u_starLR_y(bool L,ConservedStateVector& output,const double S_M,const double S_L,const double S_R,const PrimitiveStateVector& pL,const PrimitiveStateVector& pR,double gamma){
    double TOL=1e-8; //for detecting 0 denominator
    const PrimitiveStateVector& p = (L ? pL : pR);
    double S = (L ? S_L : S_R);

    output.density()=p.density()*(S-p.velocity().y())/(S-S_M);
    output.momentum().y()=output.density()*S_M;
    
    double vel_denom= (p.density()*(S-p.velocity().y())*(S-S_M)-p.B().y()*p.B().y());
    if (vel_denom<TOL){
        output.momentum().x()=output.density()*p.velocity().x();
        output.momentum().z()=output.density()*p.velocity().z();
    }
    else{
    output.momentum().x()= output.density()*(p.velocity().x()-((S_M-p.velocity().y())*p.B().y()*p.B().x())/vel_denom);
    output.momentum().z()= output.density()*(p.velocity().z()-((S_M-p.velocity().y())*p.B().y()*p.B().z())/vel_denom);
    }
    output.B().y()=p.B().y();

    double B_denom=(p.density()*(S-p.velocity().y())*(S-S_M)-p.B().y()*p.B().y());
    if (B_denom<TOL){
        output.B().x()=0.;
        output.B().z()=0.;
    }
    else{
    output.B().x()=p.B().x()*(p.density()*(S-p.velocity().y())*(S-p.velocity().y())-p.B().y()*p.B().y())/B_denom;
    output.B().z()=p.B().z()*(p.density()*(S-p.velocity().y())*(S-p.velocity().y())-p.B().y()*p.B().y())/B_denom;
    }

    double p_T_star=((S_R-pR.velocity().y())*pR.density()*pL.pressure_T()-(S_L-pL.velocity().y())*pL.density()*pR.pressure_T()+pL.density()*pR.density()*(S_R-pR.velocity().y())*(S_L-pL.velocity().y())*(pR.velocity().y()-pL.velocity().y()))/((S_R-pR.velocity().y())*pR.density()-(S_L-pL.velocity().y())*pL.density());
    //Intermediate step for convenience
    double E=p.pressure()/(gamma-1)+0.5*p.density()*dot(p.velocity(),p.velocity())+ 0.5* dot(p.B(),p.B());

    output.energy()=((S-p.velocity().y())*E-p.pressure_T()*p.velocity().y()+p_T_star*S_M + p.B().y()*(dot(p.velocity(),p.B())-1/output.density()*dot(output.momentum(),output.B())))/(S-S_M);
    output.psi()=p.psi();

}

void calc_u_doublestarLR_x(bool L,ConservedStateVector& output,ConservedStateVector& u_starL,ConservedStateVector& u_starR){
    const ConservedStateVector& u = (L ? u_starL : u_starR);
    int pm =(L? -1:1);

    output.density()=u.density();
    output.momentum().x()=u.momentum().x();
    
    int sgnBx=(0.0 < u.B().x()) - (u.B().x() < 0.0);
    
    output.momentum().y()=output.density()*(u_starL.momentum().y()/std::sqrt(u_starL.density())+u_starR.momentum().y()/std::sqrt(u_starR.density())+(u_starR.B().y()-u_starL.B().y())*sgnBx)/(std::sqrt(u_starR.density())+std::sqrt(u_starL.density()));
    output.momentum().z()=output.density()*(u_starL.momentum().z()/std::sqrt(u_starL.density())+u_starR.momentum().z()/std::sqrt(u_starR.density())+(u_starR.B().z()-u_starL.B().z())*sgnBx)/(std::sqrt(u_starR.density())+std::sqrt(u_starL.density()));

    output.B().x()=u.B().x();

    output.B().y()=(std::sqrt(u_starL.density())*u_starR.B().y()+std::sqrt(u_starR.density())*u_starL.B().y()+std::sqrt(u_starR.density()*u_starL.density())*(u_starR.momentum().y()/u_starR.density()-u_starL.momentum().y()/u_starL.density())*sgnBx)/(std::sqrt(u_starR.density())+std::sqrt(u_starL.density()));
    
    output.B().z()=(std::sqrt(u_starL.density())*u_starR.B().z()+std::sqrt(u_starR.density())*u_starL.B().z()+std::sqrt(u_starR.density()*u_starL.density())*(u_starR.momentum().z()/u_starR.density()-u_starL.momentum().z()/u_starL.density())*sgnBx)/(std::sqrt(u_starR.density())+std::sqrt(u_starL.density()));

    output.energy()=u.energy()+pm*1/std::sqrt(u.density())*(dot(u.momentum(),u.B())-dot(output.momentum(),output.B()))*sgnBx;

    output.psi()=u.psi();
}

void calc_u_doublestarLR_y(bool L,ConservedStateVector& output,ConservedStateVector& u_starL,ConservedStateVector& u_starR){
    const ConservedStateVector& u = (L ? u_starL : u_starR);
    int pm =(L? -1:1);

    output.density()=u.density();
    output.momentum().y()=u.momentum().y();
    
    int sgnBy=(0.0 < u.B().y()) - (u.B().y() < 0.0);
    
    output.momentum().x()=output.density()*(u_starL.momentum().x()/std::sqrt(u_starL.density())+u_starR.momentum().x()/std::sqrt(u_starR.density())+(u_starR.B().x()-u_starL.B().x())*sgnBy)/(std::sqrt(u_starR.density())+std::sqrt(u_starL.density()));
    output.momentum().z()=output.density()*(u_starL.momentum().z()/std::sqrt(u_starL.density())+u_starR.momentum().z()/std::sqrt(u_starR.density())+(u_starR.B().z()-u_starL.B().z())*sgnBy)/(std::sqrt(u_starR.density())+std::sqrt(u_starL.density()));

    output.B().y()=u.B().y();

    output.B().x()=(std::sqrt(u_starL.density())*u_starR.B().x()+std::sqrt(u_starR.density())*u_starL.B().x()+std::sqrt(u_starR.density()*u_starL.density())*(u_starR.momentum().x()/u_starR.density()-u_starL.momentum().x()/u_starL.density())*sgnBy)/(std::sqrt(u_starR.density())+std::sqrt(u_starL.density()));
    output.B().z()=(std::sqrt(u_starL.density())*u_starR.B().z()+std::sqrt(u_starR.density())*u_starL.B().z()+std::sqrt(u_starR.density()*u_starL.density())*(u_starR.momentum().z()/u_starR.density()-u_starL.momentum().z()/u_starL.density())*sgnBy)/(std::sqrt(u_starR.density())+std::sqrt(u_starL.density()));

    output.energy()=u.energy()+pm*1/std::sqrt(u.density())*(dot(u.momentum(),u.B())-dot(output.momentum(),output.B()))*sgnBy;

    output.psi()=u.psi();
}

void calc_HLLD_flux_x(StateVector& flux_out,const StateVector& flux_L,const StateVector& flux_R, ConservedStateVector& u_starL,ConservedStateVector& u_starR,ConservedStateVector& u_doublestar,const PrimitiveStateVector& pL,const PrimitiveStateVector& pR,const ConservedStateVector& uL,const ConservedStateVector& uR,const double gamma){
    Region region;
    auto [S_L,S_R]= calc_S_LR_x(pL,pR,gamma);
    double S_M=0;
    double S_starL=0;
    double S_starR=0;
    if(S_L>0){
        region=Region::L;
    }
    else if(S_R<0){
        region=Region::R;
    }
    else{
        S_M = calc_S_M_x(S_L,S_R,pL,pR);
        auto [S_starL,S_starR] = calc_S_stars_x(S_L,S_R,S_M,pL,pR);
        if(S_starL>0){
            region=Region::L_star;
            calc_u_starLR_x(true,u_starL,S_M,S_L,S_R,pL,pR,gamma);
        }else if(S_M>0){
            region=Region::L_doublestar;
            calc_u_starLR_x(true,u_starL,S_M,S_L,S_R,pL,pR,gamma);//left star states
            calc_u_starLR_x(false,u_starR,S_M,S_L,S_R,pL,pR,gamma);//right star states

            calc_u_doublestarLR_x(true,u_doublestar,u_starL,u_starR); //Calculate left version of double star state
        }else if(S_starR>0){
            region=Region::R_doublestar;
            calc_u_starLR_x(true,u_starL,S_M,S_L,S_R,pL,pR,gamma);//left star states
            calc_u_starLR_x(false,u_starR,S_M,S_L,S_R,pL,pR,gamma);//right star states

            calc_u_doublestarLR_x(false,u_doublestar,u_starL,u_starR);

        }else{
            region=Region::R_star;
            calc_u_starLR_x(false,u_starR,S_M,S_L,S_R,pL,pR,gamma);

        }
    }
    switch(region){
        case Region::L:
        flux_out=flux_L;
        break;
        case Region::L_star:
        flux_out+=(flux_L+S_L*(u_starL-uL));
        break;
        case Region::R_star:
        flux_out+=(flux_R+S_R*(u_starR-uR));
        break;
        case Region::L_doublestar:
        flux_out+=(flux_L+S_L*(u_starL-uL)+S_starL*(u_doublestar-u_starL));
        break;
        case Region::R_doublestar:
        flux_out+=(flux_R+S_R*(u_starR-uR)+S_starR*(u_doublestar-u_starR));
        break;
        case Region::R:
        flux_out=flux_R;
        break;

    }  
}


void calc_HLLD_flux_y(StateVector& flux_out,const StateVector& flux_L,const StateVector& flux_R, ConservedStateVector& u_starL,ConservedStateVector& u_starR,ConservedStateVector& u_doublestar,const PrimitiveStateVector& pL,const PrimitiveStateVector& pR,const ConservedStateVector& uL,const ConservedStateVector& uR,const double gamma){
    Region region;
    auto [S_L,S_R]= calc_S_LR_y(pL,pR,gamma);
    double S_M=0;
    double S_starL=0;
    double S_starR=0;
    if(S_L>0){
        region=Region::L;
    }
    else if(S_R<0){
        region=Region::R;
    }
    else{
        S_M = calc_S_M_y(S_L,S_R,pL,pR);
        auto [S_starL,S_starR] = calc_S_stars_y(S_L,S_R,S_M,pL,pR);
        if(S_starL>0){
            region=Region::L_star;
            calc_u_starLR_y(true,u_starL,S_M,S_L,S_R,pL,pR,gamma);
        }else if(S_M>0){
            region=Region::L_doublestar;
            calc_u_starLR_y(true,u_starL,S_M,S_L,S_R,pL,pR,gamma);//left star states
            calc_u_starLR_y(false,u_starR,S_M,S_L,S_R,pL,pR,gamma);//right star states

            calc_u_doublestarLR_y(true,u_doublestar,u_starL,u_starR); //Calculate left version of double star state
        }else if(S_starR>0){
            region=Region::R_doublestar;
            calc_u_starLR_y(true,u_starL,S_M,S_L,S_R,pL,pR,gamma);//left star states
            calc_u_starLR_y(false,u_starR,S_M,S_L,S_R,pL,pR,gamma);//right star states

            calc_u_doublestarLR_y(false,u_doublestar,u_starL,u_starR);

        }else{
            region=Region::R_star;
            calc_u_starLR_y(false,u_starR,S_M,S_L,S_R,pL,pR,gamma);

        }
    }
    switch(region){
        case Region::L:
        flux_out+=flux_L;
        break;
        case Region::L_star:
        flux_out+=(flux_L+S_L*(u_starL-uL));
        break;
        case Region::R_star:
        flux_out+=(flux_R+S_R*(u_starR-uR));
        break;
        case Region::L_doublestar:
        flux_out+=(flux_L+S_L*(u_starL-uL)+S_starL*(u_doublestar-u_starL));
        break;
        case Region::R_doublestar:
        flux_out+=(flux_R+S_R*(u_starR-uR)+S_starR*(u_doublestar-u_starR));
        break;
        case Region::R:
        flux_out+=flux_R;
        break;
    }  
}

void do_HLLD_x_update(Grid& grid, const SimulationConfig& cfg){
    //Starting with ubarplus states after SLIC method
    redefine_LRx(grid,cfg);
    set_xtilde_vals(grid);
    set_u_bar(grid,cfg);
    update_fluxes_x(grid);

    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    for(size_t i=1;i<nx+2*g-2;i++){
        for(size_t j=1;j<ny+2*g-2;j++){
            calc_HLLD_flux_x(grid.fluxHLLD(i,j),grid.fluxL(i,j),grid.fluxR(i,j),grid.uStarL(i,j),grid.uStarR(i,j),grid.uDoubleStar(i,j),grid.primL(i,j),grid.primR(i,j),grid.uBarL(i,j),grid.uBarR(i,j),cfg.gamma);

}}
 //update all real cells   
for(size_t i=2;i<nx+2*g-2;i++){
        for(size_t j=2;j<ny+2*g-2;j++){
            grid.U(i,j)-=grid.dt/grid.dx*(grid.fluxHLLD(i,j)-grid.fluxHLLD(i-1,j));
        }}
update_bcs(grid,cfg,grid.U);
for(size_t i=0;i<nx+2*g;i++){
        for(size_t j=0;j<ny+2*g;j++){
            grid.Prim(i,j)=grid.U(i,j).con_to_prim(cfg.gamma);//Updates both primitive and conserved data
}}}

void do_HLLD_y_update(Grid& grid, const SimulationConfig& cfg){
    //Starting with ubarplus states after SLIC method
    redefine_LRy(grid,cfg);
    set_ytilde_vals(grid);
    set_u_bar(grid,cfg);
    update_fluxes_y(grid);

    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;
    for(size_t i=1;i<nx+2*g-2;i++){
        for(size_t j=1;j<ny+2*g-2;j++){
            calc_HLLD_flux_y(grid.fluxHLLD(i,j),grid.fluxL(i,j),grid.fluxR(i,j),grid.uStarL(i,j),grid.uStarR(i,j),grid.uDoubleStar(i,j),grid.primL(i,j),grid.primR(i,j),grid.uBarL(i,j),grid.uBarR(i,j),cfg.gamma);

}}
 //update all real cells   
for(size_t i=2;i<nx+2*g-2;i++){
        for(size_t j=2;j<ny+2*g-2;j++){
            grid.U(i,j)-=grid.dt/grid.dx*(grid.fluxHLLD(i,j)-grid.fluxHLLD(i,j-1));
        }}
update_bcs(grid,cfg,grid.U);
for(size_t i=0;i<nx+2*g;i++){
        for(size_t j=0;j<ny+2*g;j++){
            grid.Prim(i,j)=grid.U(i,j).con_to_prim(cfg.gamma);//Updates both primitive and conserved data

}}
}