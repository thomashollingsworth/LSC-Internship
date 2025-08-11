#include "Initialise.h"
#include "BoundaryConditions.h"


void initialise(Grid& grid, SimulationConfig& cfg){
    using test=SimulationConfig::InitialCondition;

    size_t nx=grid.num_xcells;
    size_t ny=grid.num_ycells;
    size_t g=grid.ghost_cells;

    //Declaring parameters
    double x0;
    double xf;
    double y0;
    double yf;
    double tf;
    double gamma;
    double dx;
    double dy;

    //constants
    double pi = 4*std::atan(1);


    switch(cfg.initialcondition){
    //______________________________________________________________________
    case test::BrioWu_x:
        //Boundary Conditions
        cfg.bcs_x0=SimulationConfig::BoundaryCondition::Periodic;
        cfg.bcs_xf=SimulationConfig::BoundaryCondition::Periodic;
        cfg.bcs_y0=SimulationConfig::BoundaryCondition::Periodic;
        cfg.bcs_yf=SimulationConfig::BoundaryCondition::Periodic;
        //Domain conditions
         x0=0.;
         xf=800.;
         y0=0.;
         yf=800.;
         tf=80.;
        
        //Other parameters
         gamma=2.;

        //Setting values
        cfg.t_end=tf;
        cfg.gamma=gamma;
        
        //Initialising the x and y arrays

         dx=(xf-x0)/nx;
        grid.dx=dx;
        for(size_t i=0;i<nx;i++){
            grid.x[i]=x0+0.5*(1+i)*dx;
        }

         dy=(yf-y0)/nx;
        grid.dy=dy;
        for(size_t i=0;i<ny;i++){
            grid.y[i]=y0+0.5*(1+i)*dy;
        }

        //Move from large x to small x with field
        for(size_t i=g; i<nx+g;i++){
        for(size_t j=g;j<ny+g;j++){
        if(grid.x[i]<=400 && grid.y[j]<=400){
            //bottom left corner
            grid.Prim(i,j)[0]=0.125;
            grid.Prim(i,j)[1]=0.;
            grid.Prim(i,j)[2]=0.;
            grid.Prim(i,j)[3]=0.;
            grid.Prim(i,j)[4]=0.1;
            grid.Prim(i,j)[5]=0.75;
            grid.Prim(i,j)[6]=-1.;
            grid.Prim(i,j)[7]=0.0;
            grid.Prim(i,j)[8]=0.;
        }
        else if(grid.x[i]<=400 && grid.y[j]>400){
            //top left corner
            grid.Prim(i,j)[0]=0.125;
            grid.Prim(i,j)[1]=0.;
            grid.Prim(i,j)[2]=0.;
            grid.Prim(i,j)[3]=0.;
            grid.Prim(i,j)[4]=0.1;
            grid.Prim(i,j)[5]=0.75;
            grid.Prim(i,j)[6]=-1.;
            grid.Prim(i,j)[7]=0.0;
            grid.Prim(i,j)[8]=0.;
        }
        else if(grid.x[i]>400 && grid.y[j]<400){
            //bottom right corner
            grid.Prim(i,j)[0]=1.;
            grid.Prim(i,j)[1]=0.;
            grid.Prim(i,j)[2]=0.;
            grid.Prim(i,j)[3]=0.;
            grid.Prim(i,j)[4]=1.;
            grid.Prim(i,j)[5]=0.75;
            grid.Prim(i,j)[6]=1.;
            grid.Prim(i,j)[7]=0.;
            grid.Prim(i,j)[8]=0.;
            
        }else if(grid.x[i]>400 && grid.y[j]>400){
            //top right corner
            grid.Prim(i,j)[0]=1.;
            grid.Prim(i,j)[1]=0.;
            grid.Prim(i,j)[2]=0.;
            grid.Prim(i,j)[3]=0.;
            grid.Prim(i,j)[4]=1.;
            grid.Prim(i,j)[5]=0.75;
            grid.Prim(i,j)[6]=1.;
            grid.Prim(i,j)[7]=0.;
            grid.Prim(i,j)[8]=0.;
        }
        grid.U(i,j)=grid.Prim(i,j).prim_to_con(cfg.gamma);}}
        update_bcs(grid,cfg,grid.U);
        update_bcs(grid,cfg,grid.Prim);
        break;
    //_____________________________________________________________________________________
    
    case test::BrioWu_y:
        //Boundary Conditions
        cfg.bcs_x0=SimulationConfig::BoundaryCondition::Transimssive;
        cfg.bcs_xf=SimulationConfig::BoundaryCondition::Transimssive;
        cfg.bcs_y0=SimulationConfig::BoundaryCondition::Transimssive;
        cfg.bcs_yf=SimulationConfig::BoundaryCondition::Transimssive;
        //Domain conditions
         x0=0.;
         xf=800.;
         y0=0.;
         yf=800.;
         tf=80.;
        
        //Other parameters
         gamma=2.;

        //Setting values
        cfg.t_end=tf;
        cfg.gamma=gamma;
        
        //Initialising the x and y arrays

         dx=(xf-x0)/nx;
        grid.dx=dx;
        for(size_t i=0;i<nx;i++){
            grid.x[i]=x0+0.5*(1+i)*dx;
        }

         dy=(yf-y0)/nx;
        grid.dy=dy;
        for(size_t i=0;i<ny;i++){
            grid.y[i]=y0+0.5*(1+i)*dy;
        }

        //Move from large x to small x with field
        for(size_t i=g; i<nx+g;i++){
        for(size_t j=g;j<ny+g;j++){
        //Move from small y to large y with B field
        if(grid.x[i]<=400 && grid.y[j]<=400){
            //bottom left corner
            grid.Prim(i,j)[0]=1.;
            grid.Prim(i,j)[1]=0.;
            grid.Prim(i,j)[2]=0.;
            grid.Prim(i,j)[3]=0.;
            grid.Prim(i,j)[4]=1.;
            grid.Prim(i,j)[6]=0.75;
            grid.Prim(i,j)[5]=1.;
            grid.Prim(i,j)[7]=0.;
            grid.Prim(i,j)[8]=0.;
        }
        else if(grid.x[i]<=400 && grid.y[j]>400){
            //top left corner
            grid.Prim(i,j)[0]=0.125;
            grid.Prim(i,j)[1]=0.;
            grid.Prim(i,j)[2]=0.;
            grid.Prim(i,j)[3]=0.;
            grid.Prim(i,j)[4]=0.1;
            grid.Prim(i,j)[6]=0.75;
            grid.Prim(i,j)[5]=-1.;
            grid.Prim(i,j)[7]=0.0;
            grid.Prim(i,j)[8]=0.;
        }
        else if(grid.x[i]>400 && grid.y[j]<400){
            //bottom right corner
            grid.Prim(i,j)[0]=1.;
            grid.Prim(i,j)[1]=0.;
            grid.Prim(i,j)[2]=0.;
            grid.Prim(i,j)[3]=0.;
            grid.Prim(i,j)[4]=1.;
            grid.Prim(i,j)[6]=0.75;
            grid.Prim(i,j)[5]=1.;
            grid.Prim(i,j)[7]=0.;
            grid.Prim(i,j)[8]=0.;
        }else if(grid.x[i]>400 && grid.y[j]>400){
            //top right corner
            grid.Prim(i,j)[0]=0.125;
            grid.Prim(i,j)[1]=0.;
            grid.Prim(i,j)[2]=0.;
            grid.Prim(i,j)[3]=0.;
            grid.Prim(i,j)[4]=0.1;
            grid.Prim(i,j)[6]=0.75;
            grid.Prim(i,j)[5]=-1.;
            grid.Prim(i,j)[7]=0.0;
            grid.Prim(i,j)[8]=0.;
        }
        grid.U(i,j)=grid.Prim(i,j).prim_to_con(cfg.gamma);}}
        update_bcs(grid,cfg,grid.U);
        update_bcs(grid,cfg,grid.Prim);
        break;
    //_____________________________________________________________________________________
    
    case test::OrzsagTang:
        //Boundary Conditions
        cfg.bcs_x0=SimulationConfig::BoundaryCondition::Periodic;
        cfg.bcs_xf=SimulationConfig::BoundaryCondition::Periodic;
        cfg.bcs_y0=SimulationConfig::BoundaryCondition::Periodic;
        cfg.bcs_yf=SimulationConfig::BoundaryCondition::Periodic;
        //Domain conditions
         x0=0.;
         xf=1.;
         y0=0.;
         yf=1.;
         tf=1.1;
        
        //Other parameters
         gamma=5./3.;

        //Setting values
        cfg.t_end=tf;
        cfg.gamma=gamma;
        
        //Initialising the x and y arrays

         dx=(xf-x0)/nx;
        grid.dx=dx;
        for(size_t i=0;i<nx;i++){
            grid.x[i]=x0+0.5*(1+i)*dx;
        }

         dy=(yf-y0)/nx;
        grid.dy=dy;
        for(size_t i=0;i<ny;i++){
            grid.y[i]=y0+0.5*(1+i)*dy;
        }

        //Move from large x to small x with field
        
        
        for(size_t i=g; i<nx+g;i++){
        for(size_t j=g;j<ny+g;j++){
            grid.Prim(i,j)[0]=cfg.gamma*cfg.gamma;
            grid.Prim(i,j)[1]=-std::sin(2*pi*grid.y[j]);
            grid.Prim(i,j)[2]=std::sin(2*pi*grid.x[i]);
            grid.Prim(i,j)[3]=0.;
            grid.Prim(i,j)[4]=cfg.gamma;
            grid.Prim(i,j)[5]=-std::sin(2*pi*grid.y[j]);
            grid.Prim(i,j)[6]=std::sin(4*pi*grid.x[i]);
            grid.Prim(i,j)[7]=0;
            grid.Prim(i,j)[8]=0;
   
        grid.U(i,j)=grid.Prim(i,j).prim_to_con(cfg.gamma);}}
        
        update_bcs(grid,cfg,grid.U);
        update_bcs(grid,cfg,grid.Prim);
        break;
    //_____________________________________________________________________________________
    
    case test::KelvinHelmholtz:
        //Boundary Conditions
        cfg.bcs_x0=SimulationConfig::BoundaryCondition::Periodic;
        cfg.bcs_xf=SimulationConfig::BoundaryCondition::Periodic;
        cfg.bcs_y0=SimulationConfig::BoundaryCondition::Reflective;
        cfg.bcs_yf=SimulationConfig::BoundaryCondition::Reflective;
        //Domain conditions
         x0=0.;
         xf=1.;
         y0=-1.;
         yf=1.;
         tf=20.;
        
        //Other parameters
         gamma=5./3.;

        //Setting values
        cfg.t_end=tf;
        cfg.gamma=gamma;
        
        //Initialising the x and y arrays

         dx=(xf-x0)/nx;
        grid.dx=dx;
        for(size_t i=0;i<nx;i++){
            grid.x[i]=x0+0.5*(1+i)*dx;
        }

         dy=(yf-y0)/nx;
        grid.dy=dy;
        for(size_t i=0;i<ny;i++){
            grid.y[i]=y0+0.5*(1+i)*dy;
        }

        //Move from large x to small x with field
         pi = 4*std::atan(1);
        
        for(size_t i=g; i<nx+g;i++){
        for(size_t j=g;j<ny+g;j++){
            grid.Prim(i,j)[0]=cfg.gamma*cfg.gamma;
            grid.Prim(i,j)[1]=-std::sin(2*pi*grid.y[j]);
            grid.Prim(i,j)[2]=std::sin(2*pi*grid.x[i]);
            grid.Prim(i,j)[3]=0.;
            grid.Prim(i,j)[4]=cfg.gamma;
            grid.Prim(i,j)[5]=-std::sin(2*pi*grid.y[j]);
            grid.Prim(i,j)[6]=std::sin(4*pi*grid.x[i]);
            grid.Prim(i,j)[7]=0;
            grid.Prim(i,j)[8]=0;
        grid.U(i,j)=grid.Prim(i,j).prim_to_con(cfg.gamma);}}
        
        update_bcs(grid,cfg,grid.U);
        update_bcs(grid,cfg,grid.Prim);
        break;
}
}