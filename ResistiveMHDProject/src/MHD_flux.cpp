//Outputs the flux
#include "MHD_flux.h"

void MHD_xflux(StateVector& out, const ConservedStateVector& CSV, const PrimitiveStateVector& PSV,const double c_h){

    out[0]=CSV.momentum().x();
    out[1]=CSV.density()*(PSV.velocity().x()*PSV.velocity().x()) + PSV.pressure_T()-(CSV.B().x()*CSV.B().x());
    out[2]=CSV.density()*(PSV.velocity().x()*PSV.velocity().y())-CSV.B().x()*CSV.B().y();
    out[3]=CSV.density()*(PSV.velocity().x()*PSV.velocity().z())-CSV.B().x()*CSV.B().z();
    out[4]=PSV.velocity().x()*(CSV.energy()+PSV.pressure_T())-dot(PSV.velocity(),CSV.B())*CSV.B().x();
    //out[5]=CSV.psi();
    out[6]=CSV.B().y()*PSV.velocity().x()-CSV.B().x()*PSV.velocity().y();
    out[7]=CSV.B().z()*PSV.velocity().x()-CSV.B().x()*PSV.velocity().z();
    //out[8]=CSV.B().x()*c_h*c_h;
    out[5]=0;
    out[8]=0;

    

};

void MHD_yflux(StateVector& out, const ConservedStateVector& CSV, const PrimitiveStateVector& PSV,const double c_h){

    out[0]=CSV.momentum().y();
    out[1]=CSV.density()*(PSV.velocity().y()*PSV.velocity().x())-CSV.B().y()*CSV.B().x();
    out[2]=CSV.density()*(PSV.velocity().y()*PSV.velocity().y()) + PSV.pressure() + 0.5*dot(CSV.B(),CSV.B())-(CSV.B().y()*CSV.B().y());
    out[3]=CSV.density()*(PSV.velocity().y()*PSV.velocity().z())-CSV.B().y()*CSV.B().z();
    out[4]=PSV.velocity().y()*(CSV.energy()+PSV.pressure_T())-dot(PSV.velocity(),CSV.B())*CSV.B().y();
    out[5]=CSV.B().x()*PSV.velocity().y()-CSV.B().y()*PSV.velocity().x();
    out[6]=CSV.psi();
    out[7]=CSV.B().z()*PSV.velocity().y()-CSV.B().y()*PSV.velocity().z();
    out[8]=CSV.B().y()*c_h*c_h;

};