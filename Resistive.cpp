#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <fstream>
#include <string>
#include <cstdlib>
#include <filesystem>  // C++17
#include <initializer_list>

typedef std::array< double ,9> StateVector;

//________________________________________________________________________________________________________________________________________________
//Defining a few vector calculus operations that will work on my 2D grid for a specified parameter
//Each function will take the entire 2D grid and chosen parameter as an input and output a field 
//0th order terms are defined for all cells (both ghost layers)
//First order derivatives will be defined for all real cells + first layer of ghost cells
//Second order derivatives will only be defined for real cells

//Functions are overloaded to handle vectors or references (slices from larger vectors) 
//Additional utility function get_slice() returns a 2D array of pointers from the larger 3D array

std::vector<std::vector<const double*> > get_slice(const std::vector<std::vector<StateVector > >& u, const int slice_index){
    //Used to access 2D slices of u without the ability to modify the original u values
    std::vector<std::vector<const double* > > slice(u.size(),std::vector<const double*>(u[0].size()));
    for(size_t i=1; i<u.size()-1;i++){
        for(size_t j=1;j<u[0].size()-1;j++){

            slice[i][j]= &u[i][j][slice_index];

}}
return slice;
}

std::tuple< std::vector<std::vector<double > >,std::vector<std::vector<double > > > grad2D(const std::vector<std::vector<double> >& u,const double dx,const double dy){
    //Returns {grad_x(u),grad_y(u)} for all real cells on 2D grid u (which has 2 ghost cell boundary) plus the first ghost layer 
    std::vector<std::vector<double > > grad_x(u.size(),std::vector<double>(u[0].size()));
    std::vector<std::vector<double > > grad_y(u.size(),std::vector<double>(u[0].size()));

    for(size_t i=1; i<u.size()-1;i++){
        for(size_t j=1;j<u[0].size()-1;j++){

            grad_x[i][j]=1/(2*dx)*(u[i+1][j]-u[i-1][j]);
            grad_y[i][j]=1/(2*dy)*(u[i][j+1]-u[i][j-1]);
            
}
    }
    return {grad_x,grad_y};
}
std::tuple< std::vector<std::vector<double > >,std::vector<std::vector<double > > > grad2D(const std::vector<std::vector<const double*> >& u,const double dx,const double dy){
    //Returns {grad_x(u),grad_y(u)} for all real cells on 2D grid u (which has 2 ghost cell boundary) plus the first ghost layer 
    std::vector<std::vector<double > > grad_x(u.size(),std::vector<double>(u[0].size()));
    std::vector<std::vector<double > > grad_y(u.size(),std::vector<double>(u[0].size()));

    for(size_t i=1; i<u.size()-1;i++){
        for(size_t j=1;j<u[0].size()-1;j++){

            grad_x[i][j]=1/(2*dx)*(*u[i+1][j]-*u[i-1][j]);
            grad_y[i][j]=1/(2*dy)*(*u[i][j+1]-*u[i][j-1]);
            
}
    }
    return {grad_x,grad_y};
}

std::vector<std::vector<double > > div2D(const std::vector<std::vector<double > >& F_x, const std::vector<std::vector<double > >& F_y, const double dx,const double dy){
    //Returns Divergence of a vector quantity (F_x,F_y,F_z) for all real cells on 2D grid (which has 2 ghost cell boundary) plus the first ghost layer
    std::vector<std::vector<double > > div(F_x.size(),std::vector<double>(F_x[0].size()));

    for(size_t i=1; i<F_x.size()-1;i++){
        for(size_t j=1;j<F_x[0].size()-1;j++){

            div[i][j]=1/(2*dx)*(F_x[i+1][j]-F_x[i-1][j]);
            div[i][j]+=1/(2*dy)*(F_y[i][j+1]-F_y[i][j-1]);       
}
    }
    return div;
}
std::vector<std::vector<double > > div2D(const std::vector<std::vector<const double* > >& F_x, const std::vector<std::vector<const double* > >& F_y, const double dx,const double dy){
    //Returns Divergence of a vector quantity (F_x,F_y,F_z) for all real cells on 2D grid (which has 2 ghost cell boundary) plus the first ghost layer
    std::vector<std::vector<double > > div(F_x.size(),std::vector<double>(F_x[0].size()));

    for(size_t i=1; i<F_x.size()-1;i++){
        for(size_t j=1;j<F_x[0].size()-1;j++){

            div[i][j]=1/(2*dx)*(*F_x[i+1][j]-*F_x[i-1][j]);
            div[i][j]+=1/(2*dy)*(*F_y[i][j+1]-*F_y[i][j-1]);       
}
    }
    return div;
}

std::tuple< std::vector<std::vector<double > >,std::vector<std::vector<double > > ,std::vector<std::vector<double > > > curl2D(const std::vector<std::vector<double > >& F_x, const std::vector<std::vector<double > >& F_y, const std::vector<std::vector<double > >& F_z,const double dx,const double dy){
    //Returns {grad_x,grad_y} for all real cells on 2D grid u (which has 2 ghost cell boundary) plus the first ghost layer
    std::vector<std::vector<double > > curl_x(F_x.size(),std::vector<double>(F_x[0].size()));
    std::vector<std::vector<double > > curl_y(F_x.size(),std::vector<double>(F_x[0].size()));
    std::vector<std::vector<double > > curl_z(F_x.size(),std::vector<double>(F_x[0].size()));

    for(size_t i=1; i<F_x.size()-1;i++){
        for(size_t j=1;j<F_x[0].size()-1;j++){

            curl_x[i][j]=1/(2*dy)*(F_z[i][j+1]-F_z[i][j-1]);
            curl_y[i][j]=-1/(2*dx)*(F_z[i+1][j]-F_z[i-1][j]);
            curl_z[i][j]=1/(2*dx)*(F_y[i+1][j]-F_y[i-1][j])-1/(2*dy)*(F_x[i][j+1]-F_x[i][j-1]);
            
}
    }
    return {curl_x,curl_y,curl_z};
}
std::tuple< std::vector<std::vector<double > >,std::vector<std::vector<double > > ,std::vector<std::vector<double > > > curl2D(const std::vector<std::vector<const double* > >& F_x, const std::vector<std::vector<const double* > >& F_y, const std::vector<std::vector<const double* > >& F_z,const double dx,const double dy){
    //Returns {grad_x,grad_y} for all real cells on 2D grid u (which has 2 ghost cell boundary) plus the first ghost layer
    std::vector<std::vector<double > > curl_x(F_x.size(),std::vector<double>(F_x[0].size()));
    std::vector<std::vector<double > > curl_y(F_x.size(),std::vector<double>(F_x[0].size()));
    std::vector<std::vector<double > > curl_z(F_x.size(),std::vector<double>(F_x[0].size()));

    for(size_t i=1; i<F_x.size()-1;i++){
        for(size_t j=1;j<F_x[0].size()-1;j++){

            curl_x[i][j]=1/(2*dy)*(*F_z[i][j+1]-*F_z[i][j-1]);
            curl_y[i][j]=-1/(2*dx)*(*F_z[i+1][j]-*F_z[i-1][j]);
            curl_z[i][j]=1/(2*dx)*(*F_y[i+1][j]-*F_y[i-1][j])-1/(2*dy)*(*F_x[i][j+1]-*F_x[i][j-1]);
            
}
    }
    return {curl_x,curl_y,curl_z};
}

//Vector cross product function, overloaded to deal with reference slices or vectors
std::tuple< std::vector<std::vector<double > >,std::vector<std::vector<double > > ,std::vector<std::vector<double > > > cross2D(const std::vector<std::vector<double > >& F_x, const std::vector<std::vector<double > >& F_y, const std::vector<std::vector<double > >& F_z,const std::vector<std::vector<double > >& G_x, const std::vector<std::vector<double > >& G_y, const std::vector<std::vector<double > >& G_z){
    //Returns F cross G
    std::vector<std::vector<double > > cross_x(F_x.size(),std::vector<double>(F_x[0].size()));
    std::vector<std::vector<double > > cross_y(F_x.size(),std::vector<double>(F_x[0].size()));
    std::vector<std::vector<double > > cross_z(F_x.size(),std::vector<double>(F_x[0].size()));

    for(size_t i=1; i<F_x.size()-1;i++){
        for(size_t j=1;j<F_x[0].size()-1;j++){

            cross_x[i][j]=F_y[i][j]*G_z[i][j]-F_z[i][j]*G_y[i][j];
            cross_y[i][j]=F_z[i][j]*G_x[i][j]-F_x[i][j]*G_z[i][j];
            cross_z[i][j]=F_x[i][j]*G_y[i][j]-F_y[i][j]*G_x[i][j];
            
}
    }
    return {cross_x,cross_y,cross_z};
}
std::tuple< std::vector<std::vector<double > >,std::vector<std::vector<double > > ,std::vector<std::vector<double > > > cross2D(const std::vector<std::vector<const double* > >& F_x, const std::vector<std::vector<const double* > >& F_y, const std::vector<std::vector<const double* > >& F_z,const std::vector<std::vector<double > >& G_x, const std::vector<std::vector<double > >& G_y, const std::vector<std::vector<double > >& G_z){
    //Returns F cross G
    std::vector<std::vector<double > > cross_x(F_x.size(),std::vector<double>(F_x[0].size()));
    std::vector<std::vector<double > > cross_y(F_x.size(),std::vector<double>(F_x[0].size()));
    std::vector<std::vector<double > > cross_z(F_x.size(),std::vector<double>(F_x[0].size()));

    for(size_t i=1; i<F_x.size()-1;i++){
        for(size_t j=1;j<F_x[0].size()-1;j++){

            cross_x[i][j]=*F_y[i][j]*G_z[i][j]-*F_z[i][j]*G_y[i][j];
            cross_y[i][j]=*F_z[i][j]*G_x[i][j]-*F_x[i][j]*G_z[i][j];
            cross_z[i][j]=*F_x[i][j]*G_y[i][j]-*F_y[i][j]*G_x[i][j];
            
}
    }
    return {cross_x,cross_y,cross_z};
}
std::tuple< std::vector<std::vector<double > >,std::vector<std::vector<double > > ,std::vector<std::vector<double > > > cross2D(const std::vector<std::vector<const double* > >& F_x, const std::vector<std::vector<const double* > >& F_y, const std::vector<std::vector<const double* > >& F_z,const std::vector<std::vector<const double* > >& G_x, const std::vector<std::vector<const double* > >& G_y, const std::vector<std::vector<const double* > >& G_z){
    //Returns F cross G
    std::vector<std::vector<double > > cross_x(F_x.size(),std::vector<double>(F_x[0].size()));
    std::vector<std::vector<double > > cross_y(F_x.size(),std::vector<double>(F_x[0].size()));
    std::vector<std::vector<double > > cross_z(F_x.size(),std::vector<double>(F_x[0].size()));

    for(size_t i=1; i<F_x.size()-1;i++){
        for(size_t j=1;j<F_x[0].size()-1;j++){

            cross_x[i][j]=*F_y[i][j]* *G_z[i][j]-*F_z[i][j]* *G_y[i][j];
            cross_y[i][j]=*F_z[i][j]* *G_x[i][j]-*F_x[i][j]* *G_z[i][j];
            cross_z[i][j]=*F_x[i][j]* *G_y[i][j]-*F_y[i][j]* *G_x[i][j];
            
}
    }
    return {cross_x,cross_y,cross_z};
}
std::tuple< std::vector<std::vector<double > >,std::vector<std::vector<double > > ,std::vector<std::vector<double > > > cross2D(const std::vector<std::vector<double > >& F_x, const std::vector<std::vector<double > >& F_y, const std::vector<std::vector<double > >& F_z,const std::vector<std::vector<const double* > >& G_x, const std::vector<std::vector<const double* > >& G_y, const std::vector<std::vector<const double* > >& G_z){
    //Returns F cross G
    std::vector<std::vector<double > > cross_x(F_x.size(),std::vector<double>(F_x[0].size()));
    std::vector<std::vector<double > > cross_y(F_x.size(),std::vector<double>(F_x[0].size()));
    std::vector<std::vector<double > > cross_z(F_x.size(),std::vector<double>(F_x[0].size()));

    for(size_t i=1; i<F_x.size()-1;i++){
        for(size_t j=1;j<F_x[0].size()-1;j++){

            cross_x[i][j]=F_y[i][j]* *G_z[i][j]-F_z[i][j]* *G_y[i][j];
            cross_y[i][j]=F_z[i][j]* *G_x[i][j]-F_x[i][j]* *G_z[i][j];
            cross_z[i][j]=F_x[i][j]* *G_y[i][j]-F_y[i][j]* *G_x[i][j];
            
}
    }
    return {cross_x,cross_y,cross_z};
}

//________________________________________________________________________________________________________________________________________________



StateVector primitiveToConservative(const StateVector& prim,const double gamma){
    //Prim: rho,v_x,v_y,v_z,p,B_x,B_y,B_z (B_x is const.)
    //Conserved: rho, mom_x, mom_y, mom_z, U, B_x, B_y, B_z
    StateVector conserved;
    conserved[0]= prim[0];
    conserved[1]=prim[0]*prim[1];
    conserved[2]=prim[0]*prim[2];
    conserved[3]=prim[0]*prim[3];
    conserved[4]=prim[4]/((gamma-1)) +0.5*prim[0]*(prim[1]*prim[1]+prim[2]*prim[2]+prim[3]*prim[3])+0.5*(prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7]);
    conserved[5]=prim[5];
    conserved[6]=prim[6];
    conserved[7]=prim[7];
    conserved[8]=prim[8];

    return conserved;

}

StateVector conservativeToPrimitive(const StateVector& conserved,const double gamma){
    //Prim: rho,v_x,v_y,v_z,p,B_x,B_y,B_z (B_x is const.)
    //Conserved: rho, mom_x, mom_y, mom_z, U, B_x, B_y, B_z
    StateVector prim;
    prim[0]= conserved[0];
    prim[1]=conserved[1]/conserved[0];
    prim[2]=conserved[2]/conserved[0];
    prim[3]=conserved[3]/conserved[0];
    prim[4] = (gamma-1)*(conserved[4]-0.5*(conserved[1]*conserved[1]+conserved[2]*conserved[2]+conserved[3]*conserved[3])/conserved[0]-0.5*(conserved[5]*conserved[5]+conserved[6]*conserved[6]+conserved[7]*conserved[7]));
    prim[5]=conserved[5];
    prim[6]=conserved[6];
    prim[7]=conserved[7];
    prim[8]=conserved[8];

    return prim;

}
   
void save_to_file(const std::vector<std::vector<StateVector > >& u,const std::vector<double>& x, const double num_xcells,std::vector<double> y,const double num_ycells,const double gamma,const std::string dir,std::string filename){
    //Convert to primitive variables
    std::vector<std::vector<StateVector > > output(num_xcells+4, std::vector<StateVector >(num_ycells+4));
    for(size_t i=0;i<num_xcells+4;i++){
        for(size_t j=0;j<num_ycells+4;j++){
        output[i][j]=conservativeToPrimitive(u[i][j],gamma);
    }}

    //Output the data
    std::filesystem::create_directories(dir);
    
    filename=dir+filename;

    std::ofstream outfile(filename);

    for(size_t j=2; j<num_ycells+2;j++){
        for(size_t i=2; i<num_xcells+2;i++){
            outfile << x[i] << " " << y[j] << " " << output[i][j][0] << " "<< output[i][j][1]<< " "<< output[i][j][2]<< " "<< output[i][j][3]<< " "<< output[i][j][4]<< " "<< output[i][j][5]<< " "<< output[i][j][6]<< " "<< output[i][j][7]<< " "<< output[i][j][8]<<"\n";
}
outfile<<"\n"; //Blank line between y rows for correct gnuplot formatting


}
    std::string gnuplot_cmd = "gnuplot -e \"filename='" + filename + "'\" HLLC_plot_script.plt";
    std::system(gnuplot_cmd.c_str());
    std::cout << "Plot saved to " <<filename<< ".png" << std::endl;
    //std::filesystem::remove(filename);
}


void update_bcs(std::vector<std::vector<StateVector > >& u,const double num_xcells,const double num_ycells){
    // // Periodic 2D
    // u[0]=u[num_xcells];
    // u[1]=u[num_xcells+1];
    // u[num_xcells+2]= u[2];
    // u[num_xcells+3]= u[3];
    // for(size_t i=0;i<num_xcells+4;i++){
    //     //For every row set the four ghost cells
    //     u[i][0]=u[i][num_ycells];
    //     u[i][1]=u[i][num_ycells+1];
    //     u[i][num_ycells+2]=u[i][2];
    //     u[i][num_ycells+3] = u[i][3];

    //Kelvin-Helmholtz Test (periodic left-right reflective top-bottom)
    
    u[0]=u[num_xcells];
    u[1]=u[num_xcells+1];
    u[num_xcells+2]= u[2];
    u[num_xcells+3]= u[3];

    for (int i = 0; i < num_xcells + 2; ++i) {//reflect in u_y and B_y
    u[i][0] = u[i][3];
    u[i][1] = u[i][2]; // Bottom boundary
    u[i][num_ycells + 3] = u[i][num_ycells];
    u[i][num_ycells + 2] = u[i][num_ycells + 1]; // Top boundary
    //Tangent fields are reflected (require a minus sign)
    
    u[i][0][2] = -u[i][3][2];
    u[i][0][6] = -u[i][3][6];
    u[i][1][2] = -u[i][2][2];
    u[i][1][6] = -u[i][2][6];
    u[i][num_ycells + 2][2] = -u[i][num_ycells+1][2];
    u[i][num_ycells + 2][6] = -u[i][num_ycells+1][6];
    u[i][num_ycells + 3][2] = -u[i][num_ycells][2];
    u[i][num_ycells + 3][6] = -u[i][num_ycells][6];
}


    } 


//Set the initial t=0 condition from a given xarray

std::vector<std::vector<StateVector > > set_u0(const std::vector<double> x,const std::vector<double> y,const double gamma,const double num_xcells,const double num_ycells){
    StateVector prim;
    std::vector<std::vector<StateVector > > u0(x.size(),std::vector<StateVector >(y.size()));
    double pi = 4*std::atan(1);

    //Define intiial conditions here!
    for(size_t i=0; i<x.size();i++){
        for(size_t j=0;j<y.size();j++){
        
        //Shock Tube (1D tests)

        // //Move from large x to small x with field
        // if(x[i]<=400 && y[j]<=400){
        //     //bottom left corner
        //     prim[0]=0.125;
        //     prim[1]=0.;
        //     prim[2]=0.;
        //     prim[3]=0.;
        //     prim[4]=0.1;
        //     prim[5]=0.75;
        //     prim[6]=-1.;
        //     prim[7]=0.0;
        //     prim[8]=0.;
        // }
        // else if(x[i]<=400 && y[j]>400){
        //     //top left corner
        //     prim[0]=0.125;
        //     prim[1]=0.;
        //     prim[2]=0.;
        //     prim[3]=0.;
        //     prim[4]=0.1;
        //     prim[5]=0.75;
        //     prim[6]=-1.;
        //     prim[7]=0.0;
        //     prim[8]=0.;
        // }
        // else if(x[i]>400 && y[j]<400){
        //     //bottom right corner
        //     prim[0]=1.;
        //     prim[1]=0.;
        //     prim[2]=0.;
        //     prim[3]=0.;
        //     prim[4]=1.;
        //     prim[5]=0.75;
        //     prim[6]=1.;
        //     prim[7]=0.;
        //     prim[8]=0.;
            
        // }else if(x[i]>400 && y[j]>400){
        //     //top right corner
        //     prim[0]=1.;
        //     prim[1]=0.;
        //     prim[2]=0.;
        //     prim[3]=0.;
        //     prim[4]=1.;
        //     prim[5]=0.75;
        //     prim[6]=1.;
        //     prim[7]=0.;
        //     prim[8]=0.;
        // }
        
        
        ////Move from small y to large y no B field
        // if(x[i]<=400 && y[j]<=400){
        //     //bottom left corner
        //     prim[0]=1.;
        //     prim[1]=0.;
        //     prim[2]=0.;
        //     prim[3]=0.;
        //     prim[4]=1.;
        //     prim[6]=0.;
        //     prim[5]=0.;
        //     prim[7]=0.;
        //     prim[8]=0.;
        // }
        // else if(x[i]<=400 && y[j]>400){
        //     //top left corner
        //     prim[0]=0.125;
        //     prim[1]=0.;
        //     prim[2]=0.;
        //     prim[3]=0.;
        //     prim[4]=0.1;
        //     prim[6]=0.;
        //     prim[5]=0.;
        //     prim[7]=0.0;
        //     prim[8]=0.;
        // }
        // else if(x[i]>400 && y[j]<400){
        //     //bottom right corner
        //     prim[0]=1.;
        //     prim[1]=0.;
        //     prim[2]=0.;
        //     prim[3]=0.;
        //     prim[4]=1.;
        //     prim[6]=0.;
        //     prim[5]=0.;
        //     prim[7]=0.;
        //     prim[8]=0.;
        // }else if(x[i]>400 && y[j]>400){
        //     //top right corner
        //     prim[0]=0.125;
        //     prim[1]=0.;
        //     prim[2]=0.;
        //     prim[3]=0.;
        //     prim[4]=0.1;
        //     prim[6]=0.;
        //     prim[5]=0.;
        //     prim[7]=0.0;
        //     prim[8]=0.;
        // }


        // //Move from small y to large y with B field
        // if(x[i]<=400 && y[j]<=400){
        //     //bottom left corner
        //     prim[0]=1.;
        //     prim[1]=0.;
        //     prim[2]=0.;
        //     prim[3]=0.;
        //     prim[4]=1.;
        //     prim[6]=0.75;
        //     prim[5]=1.;
        //     prim[7]=0.;
        //     prim[8]=0.;
        // }
        // else if(x[i]<=400 && y[j]>400){
        //     //top left corner
        //     prim[0]=0.125;
        //     prim[1]=0.;
        //     prim[2]=0.;
        //     prim[3]=0.;
        //     prim[4]=0.1;
        //     prim[6]=0.75;
        //     prim[5]=-1.;
        //     prim[7]=0.0;
        //     prim[8]=0.;
        // }
        // else if(x[i]>400 && y[j]<400){
        //     //bottom right corner
        //     prim[0]=1.;
        //     prim[1]=0.;
        //     prim[2]=0.;
        //     prim[3]=0.;
        //     prim[4]=1.;
        //     prim[6]=0.75;
        //     prim[5]=1.;
        //     prim[7]=0.;
        //     prim[8]=0.;
        // }else if(x[i]>400 && y[j]>400){
        //     //top right corner
        //     prim[0]=0.125;
        //     prim[1]=0.;
        //     prim[2]=0.;
        //     prim[3]=0.;
        //     prim[4]=0.1;
        //     prim[6]=0.75;
        //     prim[5]=-1.;
        //     prim[7]=0.0;
        //     prim[8]=0.;
        // }

        
        // //Orszang-Tang vortex test
        // prim[0]=gamma*gamma;
        // prim[1]=-std::sin(2*pi*y[j]);
        // prim[2]=std::sin(2*pi*x[i]);
        // prim[3]=0.;
        // prim[4]=gamma;
        // prim[5]=-std::sin(2*pi*y[j]);
        // prim[6]=std::sin(4*pi*x[i]);
        // prim[7]=0;
        // prim[8]=0;

        // //Orszang-Tang vortex test
        // prim[0]=gamma*gamma;
        // prim[1]=-std::sin(2*pi*y[j]);
        // prim[2]=std::sin(2*pi*x[i]);
        // prim[3]=0.;
        // prim[4]=gamma;
        // prim[5]=-std::sin(2*pi*y[j]);
        // prim[6]=std::sin(4*pi*x[i]);
        // prim[7]=0;
        // prim[8]=0;

        //Kelvin-Helmholtz
        prim[0] = 1.0; // Density
        prim[1] = 0.5*std::tanh(20*y[j]); // Velocity
        prim[2] = 0.01*std::sin(2.0*pi*x[i])*std::exp(-y[j]*y[j]/(0.01));
        prim[3] = 0.0;
        prim[4] = 1.0 / gamma; // pressure
        prim[5] = 0.1*std::cos(pi/3.0); // magnetic field
        prim[6] = 0.0;
        prim[7] = 0.1*std::sin(pi/3.0);
        prim[8] = 0;



        u0[i][j]=primitiveToConservative(prim,gamma);
        }}
    update_bcs(u0,num_xcells,num_ycells);
       
    return u0;
    }
std::vector<std::vector<double > > set_eta(const std::vector<double> x,const std::vector<double> y,const double num_xcells,const double num_ycells){
    //Sets resistivity eta
    //Should be done in  tandem with set_u0, allows for resistivity to be defined on the 2D grid
    std::vector<std::vector<double > > eta(x.size(),std::vector<double >(y.size()));
    
    //UNIFORM CASE
    const double uniform_eta=1e-5;
    for(size_t i=0; i<x.size();i++){
        for(size_t j=0;j<y.size();j++){
            eta[i][j]=uniform_eta;        

}}
//BCs FOR PERIODIC IN X,Y

    eta[0]=eta[num_xcells];
    eta[1]=eta[num_xcells+1];
    eta[num_xcells+2]= eta[2];
    eta[num_xcells+3]= eta[3];
    for(size_t i=0;i<num_xcells+4;i++){
        //For every row set the four ghost cells
        eta[i][0]=eta[i][num_ycells];
        eta[i][1]=eta[i][num_ycells+1];
        eta[i][num_ycells+2]=eta[i][2];
        eta[i][num_ycells+3] = eta[i][3];

return eta;
}}   

std::vector<std::vector<StateVector > > calc_source_terms(const std::vector<std::vector<StateVector > >& u,const std::vector<std::vector<double > >& eta, const double dx,const double dy,const int num_xcells,const int num_ycells){
//Can customise the source terms for every parameter in this function
//Psi isn't given a source term as it can be updated analytically
//This should be calculated every time step / subtime step

//Source term for each parameter is 0 by default

std::vector<std::vector<StateVector > > sources(num_xcells+4, std::vector<StateVector >(num_ycells+4,StateVector{0.0}));

//Get 2D views of the 3D u array 
const std::vector<std::vector<const double* > >& Bx= get_slice(u,5);
const std::vector<std::vector<const double* > >& By=get_slice(u,6);
const std::vector<std::vector<const double* > >& Bz=get_slice(u,7);


//J=curlB Remember all 1st derivative quantities are only defined on real cells+1 ghost layer!
const auto& [Jx,Jy,Jz]=curl2D(Bx,By,Bz,dx,dy);

//This is a second order derivative so is only defined on real cells i.e. outer two layers of this array are useless/undefined
const auto& [curlJ_x,curlJ_y,curlJ_z]=curl2D(Jx,Jy,Jz,dx,dy);

//These terms will be used later so are convenient to calc. now
const auto& [grad_eta_x,grad_eta_y]= grad2D(eta,dx,dy);
//z component of grad will always be zeros
const std::vector<std::vector<double> > grad_eta_z(num_xcells+4, std::vector<double >(num_ycells+4,0.0));

const auto& [BcrossJ_x,BcrossJ_y,BcrossJ_z]=cross2D(Bx,By,Bz,Jx,Jy,Jz);

const auto& [grad_etacrossJ_x,grad_etacrossJ_y,grad_etacrossJ_z]=cross2D(grad_eta_x,grad_eta_y,grad_eta_z,Jx,Jy,Jz);


for(size_t i=2;i<num_xcells+2;i++){
        for(size_t j=2;j<num_ycells+2;j++){
        
        
        //Giving the energy (4) source term: Div(B cross eta*J)
        sources[i][j][4]=(grad_eta_x[i][j]*BcrossJ_x[i][j] + grad_eta_y[i][j]*BcrossJ_y[i][j])+
                        eta[i][j]*(
                        (Jx[i][j]*Jx[i][j]+Jy[i][j]*Jy[i][j]+Jz[i][j]*Jz[i][j])
                        + (*Bx[i][j] * curlJ_x[i][j] + *By[i][j] * curlJ_y[i][j] + *Bz[i][j] * curlJ_z[i][j])
                        );
        
        //Giving the B field (5,6,7) source terms: -curl(eta*J)

        sources[i][j][5]= -grad_etacrossJ_x[i][j] -eta[i][j]*curlJ_x[i][j];
        sources[i][j][6]= -grad_etacrossJ_y[i][j] -eta[i][j]*curlJ_y[i][j];
        sources[i][j][7]= -grad_etacrossJ_z[i][j] -eta[i][j]*curlJ_z[i][j];
        
    }}

return sources;
}

StateVector RK4_update(StateVector u, const double dt_source){
    
}

StateVector MHD_xflux(const StateVector& u0,const double gamma,const double c_h){
    //Prim: rho,v_x,v_y,v_z,p,B_x,B_y,B_z (B_x is const.)
    //Conserved: rho, mom_x, mom_y, mom_z, U, B_x, B_y, B_z
    StateVector flux;
    StateVector prim=conservativeToPrimitive(u0,gamma);
    flux[0]=prim[0]*prim[1];
    flux[1]=prim[0]*prim[1]*prim[1]+prim[4]+0.5*(prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7])-prim[5]*prim[5];
    flux[2]= prim[0]*prim[1]*prim[2]-prim[5]*prim[6];
    flux[3]= prim[0]*prim[1]*prim[3]-prim[5]*prim[7];
    flux[4]=(u0[4]+prim[4]+0.5*(prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7]))*prim[1]-(prim[1]*prim[5]+prim[2]*prim[6]+prim[3]*prim[7])*prim[5];

    flux[5]=u0[8];
    flux[8]=u0[5]*c_h*c_h;
    
    
    
    flux[6]=prim[6]*prim[1]-prim[5]*prim[2];
    flux[7]=prim[7]*prim[1]-prim[5]*prim[3];

    
    return flux;

}
StateVector MHD_yflux(const StateVector& u0,const double gamma,const double c_h){
    //Prim: rho,v_x,v_y,v_z,p,B_x,B_y,B_z (B_x is const.)
    //Conserved: rho, mom_x, mom_y, mom_z, U, B_x, B_y, B_z
    StateVector flux;
    StateVector prim=conservativeToPrimitive(u0,gamma);
    flux[0]=prim[0]*prim[2];
    flux[1]= prim[0]*prim[2]*prim[1]-prim[6]*prim[5];
    flux[2]=prim[0]*prim[2]*prim[2]+prim[4]+0.5*(prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7])-prim[6]*prim[6];
    flux[3]= prim[0]*prim[2]*prim[3]-prim[6]*prim[7];
    flux[4]=(u0[4]+prim[4]+0.5*(prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7]))*prim[2]-(prim[1]*prim[5]+prim[2]*prim[6]+prim[3]*prim[7])*prim[6];
    flux[5]=-(prim[6]*prim[1]-prim[2]*prim[5]);
    flux[6]=u0[8];
    flux[7]=-(prim[6]*prim[3]-prim[2]*prim[7]);
    
    flux[8]=u0[6]*c_h*c_h;
    
    return flux;



}



//HLL METHOD
double get_cfx(const double gamma, const StateVector& u){
    
    // u is the vector of parameters at lattice site
    const StateVector prim =conservativeToPrimitive(u,gamma);

    
    double factor= (gamma*prim[4]+u[5]*u[5]+u[6]*u[6]+u[7]*u[7])/u[0];

    double c_fx=std::sqrt(0.5*(factor+std::sqrt(factor*factor-4.0*gamma*prim[4]*u[5]*u[5]/(u[0]*u[0]))));
    return c_fx;

}
double get_cfy(const double gamma, const StateVector& u){
    
    // u is the vector of parameters at lattice site
    const StateVector prim =conservativeToPrimitive(u,gamma);
    
    double factor= (gamma*prim[4]+u[5]*u[5]+u[6]*u[6]+u[7]*u[7])/u[0];
    double c_fy=std::sqrt(0.5*(factor+std::sqrt(factor*factor-4.0*gamma*prim[4]*u[6]*u[6]/(u[0]*u[0]))));
    return c_fy;
}
double get_cfz(const double gamma, const StateVector& u){
    
    // u is the vector of parameters at lattice site
    const StateVector& prim =conservativeToPrimitive(u,gamma);
    
    double factor= (gamma*prim[4]+u[5]*u[5]+u[6]*u[6]+u[7]*u[7])/u[0];
    double c_fz=std::sqrt(0.5*(factor+std::sqrt(factor*factor-4.0*gamma*prim[4]*u[7]*u[7]/(u[0]*u[0]))));
    return c_fz;
}


std::tuple< double, double > get_S_LR_x(const StateVector& u_L,const StateVector& u_R, const double gamma){
    const StateVector prim_L=conservativeToPrimitive(u_L,gamma);
    const StateVector prim_R=conservativeToPrimitive(u_R,gamma);
    double v_xL = prim_L[1];
    double v_xR = prim_R[1];
    double c_fL=get_cfx(gamma,u_L);
    double c_fR=get_cfx(gamma,u_R);
    double S_L = std::min(v_xL,v_xR)-std::max(c_fL,c_fR);
    double S_R = std::max(v_xL,v_xR)+std::max(c_fL,c_fR);


    return {S_L,S_R};

}
std::tuple< double, double > get_S_LR_y(const StateVector& u_L,const StateVector& u_R, const double gamma){
    const StateVector prim_L=conservativeToPrimitive(u_L,gamma);
    const StateVector prim_R=conservativeToPrimitive(u_R,gamma);
    double v_yL = prim_L[2];
    double v_yR = prim_R[2];
    double c_fL=get_cfy(gamma,u_L);
    double c_fR=get_cfy(gamma,u_R);
    double S_L = std::min(v_yL,v_yR)-std::max(c_fL,c_fR);
    double S_R = std::max(v_yL,v_yR)+std::max(c_fL,c_fR);



    return {S_L,S_R};

}

StateVector get_u_HLL_x(const StateVector& u_L,const StateVector& u_R, const double gamma,const double c_h){
    const auto [S_L,S_R]= get_S_LR_x(u_L,u_R,gamma);

    const StateVector f_L=MHD_xflux(u_L,gamma, c_h);
    const StateVector f_R=MHD_xflux(u_R,gamma,c_h);
    StateVector u_HLL;

    for(size_t i=0;i<9;i++){
        u_HLL[i]=1/(S_R-S_L)*(S_R*u_R[i]-S_L*u_L[i]+f_L[i]-f_R[i]);
    }

    if(S_L>=0){
        return u_L;
    }else if(S_L<0 && S_R>0){
        return u_HLL;
    } else{
        return u_R;
    }
}
StateVector get_u_HLL_y(const StateVector& u_L,const StateVector& u_R, const double gamma,const double c_h){
    const auto [S_L,S_R]= get_S_LR_y(u_L,u_R,gamma);

    const StateVector f_L=MHD_yflux(u_L,gamma, c_h);
    const StateVector f_R=MHD_yflux(u_R,gamma, c_h);
    StateVector u_HLL;

    for(size_t i=0;i<9;i++){
        u_HLL[i]=1/(S_R-S_L)*(S_R*u_R[i]-S_L*u_L[i]+f_L[i]-f_R[i]);
    }

    if(S_L>=0){
        return u_L;
    }else if(S_L<0 && S_R>0){
        return u_HLL;
    }else{
        return u_R;
    }
}

StateVector get_u_star_x(const StateVector& u_L,const StateVector& u_R, const double gamma,const bool L,const double psi_tilde, const double S_L, const double S_R,const double c_h){
    //L bool dictates whether you are calculating u_starL or u_starR 
    //Everything here should have BxL=BxR= Bxtilde apart from the S_L and S_R calcs
    
    const StateVector u_HLL=get_u_HLL_x(u_L,u_R,gamma,c_h);
    const StateVector prim_HLL=conservativeToPrimitive(u_HLL,gamma);
    
    
    const StateVector prim_L= conservativeToPrimitive(u_L,gamma);
    const StateVector prim_R= conservativeToPrimitive(u_R,gamma);
    
    double S; // Will be used as S_L for L case and S_R for R_case
    StateVector prim; // Will be used as L for L case and R for R_case
    StateVector u;

    if(L){
        S=S_L;
        prim=prim_L;
        u=u_L;

    }else{
        S=S_R;
        prim=prim_R;
        u=u_R;

    }
    //universal quantity used for both L and R
    const double q_star= (prim_R[0]*prim_R[1]*(S_R-prim_R[1])-prim_L[0]*prim_L[1]*(S_L-prim_L[1])+prim_L[4]+0.5*(prim_L[5]*prim_L[5]+prim_L[6]*prim_L[6]+prim_L[7]*prim_L[7])-prim_R[4]-0.5*(prim_R[5]*prim_R[5]+prim_R[6]*prim_R[6]+prim_R[7]*prim_R[7])-prim_L[5]*prim_L[5]+prim_R[5]*prim_R[5])/(prim_R[0]*(S_R-prim_R[1])-prim_L[0]*(S_L-prim_L[1]));
    
    StateVector u_star;
    
    u_star[5]=u_HLL[5];
    u_star[6]=u_HLL[6];
    u_star[7]=u_HLL[7];

    //Calc. the momentum and density terms
    u_star[0]=prim[0]*(S-prim[1])/(S-q_star);
    u_star[1]=u_star[0]*q_star;
    u_star[2]=u_star[0]*prim[2]-(u_star[5]*u_star[6]-prim[5]*prim[6])/(S-q_star);
    u_star[3]=u_star[0]*prim[3]-(u_star[5]*u_star[7]-prim[5]*prim[7])/(S-q_star);


    //Calc. the E_star term
    double p_star= prim[0]*(S-prim[1])*(q_star-prim[1])+prim[4]+0.5*(prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7])-prim[5]*prim[5]+u_star[5]*u_star[5];
    
    u_star[4]=u[4]*(S-prim[1])/(S-q_star)+((p_star*q_star-(prim[4]+0.5*(prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7]))*prim[1])-(u_star[5]*(prim_HLL[5]*prim_HLL[1]+prim_HLL[6]*prim_HLL[2]+prim_HLL[7]*prim_HLL[3])-prim[5]*(prim[5]*prim[1]+prim[6]*prim[2]+prim[7]*prim[3])))/(S-q_star);

    u_star[8]=psi_tilde;
    
    return u_star;
}
StateVector get_u_star_y(const StateVector& u_L,const StateVector& u_R, const double gamma,const bool L,const double psi_tilde, const double S_L,const double S_R,const double c_h){
    //L bool dictates whether you are calculating u_starL or u_starR 
    
    const StateVector u_HLL=get_u_HLL_y(u_L,u_R,gamma,c_h);
    const StateVector prim_HLL=conservativeToPrimitive(u_HLL,gamma);
   
    const StateVector prim_L= conservativeToPrimitive(u_L,gamma);
    const StateVector prim_R= conservativeToPrimitive(u_R,gamma);
    
    double S; // Will be used as S_L for L case and S_R for R_case
    StateVector prim; // Will be used as L for L case and R for R_case
    StateVector u;

    if(L){
        S=S_L;
        prim=prim_L;
        u=u_L;

    }else{
        S=S_R;
        prim=prim_R;
        u=u_R;

    }
    //universal quantity used for both L and R
    const double q_star= (prim_R[0]*prim_R[2]*(S_R-prim_R[2])-prim_L[0]*prim_L[2]*(S_L-prim_L[2])+prim_L[4]+0.5*(prim_L[5]*prim_L[5]+prim_L[6]*prim_L[6]+prim_L[7]*prim_L[7])-prim_R[4]-0.5*(prim_R[5]*prim_R[5]+prim_R[6]*prim_R[6]+prim_R[7]*prim_R[7])-prim_L[6]*prim_L[6]+prim_R[6]*prim_R[6])/(prim_R[0]*(S_R-prim_R[2])-prim_L[0]*(S_L-prim_L[2]));
    
    StateVector u_star;
    
    u_star[5]=u_HLL[5];
    u_star[6]=u_HLL[6];
    u_star[7]=u_HLL[7];

    //Calc. the momentum and density terms
    u_star[0]=prim[0]*(S-prim[2])/(S-q_star);
    u_star[1]=u_star[0]*prim[1]-(u_star[5]*u_star[6]-prim[5]*prim[6])/(S-q_star);
    u_star[2]=u_star[0]*q_star;
    u_star[3]=u_star[0]*prim[3]-(u_star[6]*u_star[7]-prim[6]*prim[7])/(S-q_star);

    //Calc. the E_star term
    double p_star= prim[0]*(S-prim[2])*(q_star-prim[2])+prim[4]+0.5*(prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7])-prim[6]*prim[6]+u_star[6]*u_star[6];
    
    u_star[4]=u[4]*(S-prim[2])/(S-q_star)+((p_star*q_star-(prim[4]+0.5*(prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7]))*prim[2])-(u_star[6]*(prim_HLL[5]*prim_HLL[1]+prim_HLL[6]*prim_HLL[2]+prim_HLL[7]*prim_HLL[3])-prim[6]*(prim[5]*prim[1]+prim[6]*prim[2]+prim[7]*prim[3])))/(S-q_star);
    u_star[8]=psi_tilde;
    return u_star;
}


double get_c_h(const std::vector<std::vector<StateVector > >& u,double num_xcells,double num_ycells,double gamma){

    double max_ch=0.;
    StateVector prim;
    
    for(size_t i=2; i<num_xcells+2;i++){
        for(size_t j=2; j<num_ycells+2;j++){
        prim=conservativeToPrimitive(u[i][j],gamma);

        double c_fx=get_cfx(gamma,u[i][j]);
        double c_fy = get_cfy(gamma,u[i][j]);
        double c_fz = get_cfz(gamma,u[i][j]);

        double new_max=std::max({std::abs(prim[1])+c_fx,std::abs(prim[2])+c_fy,std::abs(prim[3])+c_fz});

        if(new_max> max_ch){
            max_ch=new_max;
        }   
    }
}

return max_ch;
}

double computeTimeStep(const double dx,const double dy, const double C,const double c_h){
 //Making an adjustment here
return C* std::min(dx,dy)/(c_h);
}


std::tuple< double,double > get_xtilde_vals(const StateVector& uL, const StateVector& uR,const double gamma,const double c_h){
    
    double Bx_tilde = 0.5*(uL[5]+uR[5])-0.5*1./c_h*(uR[8]-uL[8]);
    double psi_tilde= 0.5*(uL[8]+uR[8])-c_h/2.*(uR[5]-uL[5]);
    return {Bx_tilde,psi_tilde};

}

std::tuple< double,double > get_ytilde_vals(const StateVector& uL, const StateVector& uR,const double gamma,const double c_h){
    
    double By_tilde = 0.5*(uL[6]+uR[6])-0.5*1./c_h*(uR[8]-uL[8]);
    double psi_tilde= 0.5*(uL[8]+uR[8])-c_h/2.*(uR[6]-uL[6]);
    return {By_tilde,psi_tilde};

}

StateVector get_HLLC_flux_x(StateVector u_L,StateVector u_R, const double gamma,const double c_h){
    //gets flux between L and R position
    
    
    const auto[B_tilde,psi_tilde]=get_xtilde_vals(u_L,u_R,gamma,c_h);
    const auto [S_L,S_R]= get_S_LR_x(u_L,u_R,gamma);

    u_L[5]=B_tilde;
    u_R[5]=B_tilde;
    u_L[8] = psi_tilde;
	u_R[8] = psi_tilde;
    
    
    const StateVector f_L=MHD_xflux(u_L,gamma,c_h);
    const StateVector f_R=MHD_xflux(u_R,gamma,c_h);
    
    const StateVector u_starL=get_u_star_x(u_L,u_R,gamma,true,psi_tilde, S_L,S_R,c_h);
    const StateVector u_starR=get_u_star_x(u_L,u_R,gamma,false,psi_tilde,S_L,S_R,c_h);

    const double q_star= u_starL[1]/u_starL[0];
    
    StateVector flux_out;

    //Deciding what region/flux is being used
    if(S_L>=0){
        flux_out=f_L;

    }else if(S_L<0 && q_star>=0){
        for(size_t i=0;i<9;i++){
        flux_out[i]=f_L[i]+S_L*(u_starL[i]-u_L[i]);
    }}else if(q_star<=0 && 0<=S_R){
        for(size_t i=0;i<9;i++){
        flux_out[i]=f_R[i]+S_R*(u_starR[i]-u_R[i]);

    }}else {
        flux_out=f_R; 
    }
    
    

    // flux_out[5]=0;
    // flux_out[8]=0;

    

    return flux_out;


    // std::vector<double> fluxHLL(9);
    // for(size_t i=0;i<9;i++){
    //     fluxHLL[i]=(S_R*f_L[i]-S_L*f_R[i] +S_L*S_R*(u_R[i]-u_L[i]))/(S_R-S_L);
    // }
    // fluxHLL[5]=psi_tilde;
    // fluxHLL[8]=B_tilde*c_h*c_h;
    // return fluxHLL;
}
StateVector get_HLLC_flux_y(StateVector u_L,StateVector u_R, const double gamma,const double c_h){
    //gets flux between L and R position
    const auto[B_tilde,psi_tilde]=get_ytilde_vals(u_L,u_R,gamma,c_h);
    const auto [S_L,S_R]= get_S_LR_y(u_L,u_R,gamma);

    u_L[6]=B_tilde;
    u_R[6]=B_tilde;
    u_L[8] = psi_tilde;
	u_R[8] = psi_tilde;
    
    const StateVector f_L=MHD_yflux(u_L,gamma,c_h);
    const StateVector f_R=MHD_yflux(u_R,gamma,c_h);
    const StateVector u_starL=get_u_star_y(u_L,u_R,gamma,true,psi_tilde,S_L,S_R,c_h);
    const StateVector u_starR=get_u_star_y(u_L,u_R,gamma,false,psi_tilde, S_L,S_R,c_h);

    const double q_star= u_starL[2]/u_starL[0];
    StateVector flux_out;

    //Deciding what region/flux is being used
    if(S_L>=0){
        flux_out=f_L;
    }else if(S_L<0 && q_star>=0){
        for(size_t i=0;i<9;i++){
        flux_out[i]=f_L[i]+S_L*(u_starL[i]-u_L[i]);
    }}else if(q_star<=0 && 0<=S_R){
        for(size_t i=0;i<9;i++){
        flux_out[i]=f_R[i]+S_R*(u_starR[i]-u_R[i]);
    }}else {
        flux_out=f_R; 
    }


    flux_out[6]=psi_tilde;
    flux_out[8]=B_tilde*c_h*c_h;

    // flux_out[6]=0;
    // flux_out[8]=0;

    return flux_out;

    // std::vector<double> fluxHLL(9);
    // for(size_t i=0;i<9;i++){
    //     fluxHLL[i]=(S_R*f_L[i]-S_L*f_R[i] +S_L*S_R*(u_R[i]-u_L[i]))/(S_R-S_L);
    // }
    // fluxHLL[6]=psi_tilde;
    // fluxHLL[8]=B_tilde*c_h*c_h;
    // return fluxHLL;
}


//Adding Slope Limiting

StateVector getDelta_plus(const StateVector& u_0, const StateVector& u_plus){

    StateVector delta_plus;
    for(size_t i=0;i<9;i++){
        delta_plus[i]=u_plus[i]-u_0[i];
    }
    return delta_plus;
}
StateVector getDelta_minus(const StateVector& u_minus, const StateVector& u_0){
    StateVector delta_minus;
    for(size_t i=0;i<9;i++){
        delta_minus[i]=u_0[i]-u_minus[i];
    }
    return delta_minus;
}
StateVector getDelta(const StateVector& u_minus,const StateVector& u_0, const StateVector& u_plus,const double w=0){
StateVector delta_minus = getDelta_minus(u_minus,u_0);
StateVector delta_plus = getDelta_plus(u_0,u_plus);
StateVector delta;
    for(size_t i=0;i<9;i++){
        delta[i]=0.5*(1+w)*delta_minus[i] + 0.5*(1-w)*delta_plus[i];
    }
    return delta;
}

StateVector get_r(const StateVector& u_minus,const StateVector& u_0, const StateVector& u_plus,const double gamma){
    const StateVector delta_minus = getDelta_minus(u_minus,u_0);
    const StateVector delta_plus = getDelta_plus(u_0,u_plus);
    StateVector r;
    
    
    double min_r=2;   
    
    for(size_t i=0;i<9;i++){
    
        r[i]=delta_minus[i]/(delta_plus[i]+1e-12);
        if(r[i]<min_r && r[i] != 0.){
            min_r=r[i];
        }
        
    }
    
    // for(size_t i=0;i<9;i++){
    //     r[i]=min_r;
    //     }
        return r;
}
StateVector getXi_VL(const StateVector& r){
    //Van Leer slope limiting
    StateVector xi;

    for(size_t i=0;i<9;i++){

        if(r[i]<=0){
            xi[i]=0;
        } else{
            xi[i]=std::min(2*r[i]/(1+r[i]),2/(1+r[i]));
        } 
    }
    return xi;
}

std::tuple< StateVector,StateVector  > calc_ubar(const StateVector& u_minus,const StateVector& u_0, const StateVector& u_plus,const double w,const double gamma){
    //returns {ubar_L,ubar_R}
    const StateVector prim_minus = conservativeToPrimitive(u_minus,gamma);
    const StateVector prim_0 = conservativeToPrimitive(u_0,gamma);
    const StateVector prim_plus = conservativeToPrimitive(u_plus,gamma);
    
    StateVector ubar_L;
    StateVector ubar_R;

    const StateVector delta=getDelta(u_minus,u_0,u_plus,w);
    const StateVector r= get_r(prim_minus,prim_0,prim_plus,gamma);
    const StateVector xi= getXi_VL(r);
    
    for(size_t j=0;j<9;j++){
        //Looping over variables
            ubar_L[j]=u_0[j]-0.5*xi[j]*delta[j];
            ubar_R[j]=u_0[j]+0.5*xi[j]*delta[j];
        }  

        return {ubar_L,ubar_R};
    }
std::tuple< StateVector,StateVector > calc_ubarx_plus(const StateVector& u_minus,const StateVector& u_0, const StateVector& u_plus,const double w,const double dt, const double dx,const double gamma,const double c_h){
    //returns {ubar_L^(n+1/2),ubar_R^(n+1/2) }
    const auto [ubar_L, ubar_R]= calc_ubar(u_minus,u_0,u_plus,w,gamma);
    StateVector ubar_L_plus;
    StateVector ubar_R_plus;
    const StateVector f_R= MHD_xflux(ubar_R,gamma,c_h);
    const StateVector f_L= MHD_xflux(ubar_L,gamma,c_h);
    for(size_t j=0;j<9;j++){//Looping over variables
            ubar_L_plus[j]= ubar_L[j]-0.5*dt/dx*(f_R[j]-f_L[j]);
            ubar_R_plus[j]= ubar_R[j]-0.5*dt/dx*(f_R[j]-f_L[j]);
            
        }  
        return {ubar_L_plus,ubar_R_plus};

    }
std::tuple< StateVector,StateVector > calc_ubary_plus(const StateVector& u_minus,const StateVector& u_0, const StateVector& u_plus,const double w,const double dt, const double dx,const double gamma,const double c_h){
    //returns {ubar_L^(n+1/2),ubar_R^(n+1/2) }
    const auto [ubar_L, ubar_R]= calc_ubar(u_minus,u_0,u_plus,w,gamma);
    StateVector ubar_L_plus;
    StateVector ubar_R_plus;
    const StateVector f_R= MHD_yflux(ubar_R,gamma,c_h);
    const StateVector f_L= MHD_yflux(ubar_L,gamma,c_h);
    for(size_t j=0;j<9;j++){//Looping over variables
            ubar_L_plus[j]= ubar_L[j]-0.5*dt/dx*(f_R[j]-f_L[j]);
            ubar_R_plus[j]= ubar_R[j]-0.5*dt/dx*(f_R[j]-f_L[j]);
            
        }  
        return {ubar_L_plus,ubar_R_plus};

    }

std::tuple< std::vector<std::vector<StateVector > >,std::vector<std::vector<StateVector > > > get_ubar_plus_xarrays(const std::vector<std::vector<StateVector > >& u, const double w,const double dt, const double dx,const double gamma,const double num_xcells,const double num_ycells,const double c_h){
    std::vector<std::vector<StateVector > > ubar_Rplus(num_xcells+4, std::vector<StateVector >(num_ycells+4));
    std::vector<std::vector<StateVector > > ubar_Lplus(num_xcells+4, std::vector<StateVector >(num_ycells+4));
    
    for(size_t i=1;i<num_xcells+3;i++){
        for(size_t j=1;j<num_ycells+3;j++){
        const auto [L, R] = calc_ubarx_plus(u[i-1][j],u[i][j],u[i+1][j],w,dt,dx,gamma,c_h);
        ubar_Lplus[i][j]=L;
        ubar_Rplus[i][j]=R;  
    }}


    update_bcs(ubar_Lplus,num_xcells,num_ycells);
    update_bcs(ubar_Rplus,num_xcells,num_ycells);
    
    return {ubar_Lplus,ubar_Rplus};
}
std::tuple< std::vector<std::vector<StateVector > >,std::vector<std::vector<StateVector > > > get_ubar_plus_yarrays(const std::vector<std::vector<StateVector > >& u, const double w,const double dt, const double dx,const double gamma,const double num_xcells,const double num_ycells,const double c_h){
    std::vector<std::vector<StateVector > > ubar_Rplus(num_xcells+4, std::vector<StateVector >(num_ycells+4));
    std::vector<std::vector<StateVector > > ubar_Lplus(num_xcells+4, std::vector<StateVector >(num_ycells+4));
    
    for(size_t i=1;i<num_xcells+3;i++){
        for(size_t j=1;j<num_ycells+3;j++){
        const auto [L, R] = calc_ubary_plus(u[i][j-1],u[i][j],u[i][j+1],w,dt,dx,gamma,c_h);
        ubar_Lplus[i][j]=L;
        ubar_Rplus[i][j]=R;  
    }}

    
    update_bcs(ubar_Lplus,num_xcells,num_ycells);
    update_bcs(ubar_Rplus,num_xcells,num_ycells);
    
    return {ubar_Lplus,ubar_Rplus};
}


int main(void){
    
    // //Orszang-Tang vortex test
    // const double C=0.8;
    // const double tEnd=1.1;
    // const double tStart=0;
    // const double x0=0;
    // const double xf=1;
    // const double y0=0;
    // const double yf=1;
    // const double num_xcells=256;
    // const double num_ycells=256;
    // const double gamma=5./3. ;
    // const double w=0;
    // const int save_interval=10;

    // // Shock Tube Test
    // const double C=0.8;
    // const double tEnd=80;
    // const double tStart=0;
    // const double x0=0;
    // const double xf=800;
    // const double y0=0;
    // const double yf=800;
    // const double num_xcells=250;
    // const double num_ycells=250;
    // const double gamma=2. ;
    // const double w=0;
    // const int save_interval=10;


    //Kelvin-Helmholtz test
    const double C=0.8;
    const double tEnd=20;
    const double tStart=0;
    const double x0=0;
    const double xf=1;
    const double y0=-1;
    const double yf=1;
    const double num_xcells=128;
    const double num_ycells=256;
    const double gamma=5./3. ;
    const double w=0;
    const int save_interval=50;


    


    const double dx=(xf-x0)/num_xcells;
    const double dy=(yf-y0)/num_ycells;

    //Initialise x and y vectors which mark positions of cell centres
    std::vector<double> x(num_xcells+4);
    for(int i=0;i<x.size();i++){
        x[i]=x0+(i-0.5)*dx;

    };
    
    std::vector<double> y(num_ycells+4);
    for(int i=0;i<y.size();i++){
        y[i]=y0+(i-0.5)*dy;

    };



    //Initialise u(x) vector
    std::vector<std::vector<StateVector > > u(num_xcells+4, std::vector<StateVector >(num_ycells+4));

    // store flux: flux[i] is flux at x(i+1/2)
    std::vector<std::vector<StateVector > > flux_x(num_xcells+4, std::vector<StateVector >(num_ycells+4));
    std::vector<std::vector<StateVector > > flux_y(num_xcells+4, std::vector<StateVector >(num_ycells+4));

    u=set_u0(x,y,gamma,num_xcells,num_ycells);

    double t=tStart;
    
    int count=0; //counts number of update steps performed

    do {

        if(count%save_interval==0){
            std::string name= "Plot_" + std::to_string(t);
             //Output the data
            std::string dir = "KelvinHelmholtz/";
            save_to_file(u,x,num_xcells,y,num_ycells,gamma,dir,name);
        }

        update_bcs(u,num_xcells,num_ycells);
        
        double c_h= get_c_h(u,num_xcells,num_ycells,gamma); //Max across entirety of u, should be updated every time u is actually updated 
        double c_p_squared=c_h*0.18;
        double dt=computeTimeStep(dx,dy,C,c_h);

        t+=dt;
        count++;
        
        std::cout<<"Starting iteration: "<<count<<std::endl;

        //Do x updates first
        
        //Do subcycling time step source update for all variables with source terms updating real cells only
        
        for(size_t i=2; i<num_xcells+2; i++){
            for(size_t j=2; j<num_ycells+2; j++){

                u[i][j][8]*=std::exp((-dt/2)*(c_h*c_h)/c_p_squared);//Has exact solution so can be solved analytically at each step

                }}
        
        //First update the u_bar arrays to correspond to (x) flux for slope limiting
        const auto& [ubar_L_plusx,ubar_R_plusx] = get_ubar_plus_xarrays(u,w,dt,dx,gamma,num_xcells,num_ycells,c_h);

         // update x flux array
        

        for(size_t i=1; i<num_xcells+2; i++){
            for(size_t j=1; j<num_ycells+2; j++){

                flux_x[i][j]=get_HLLC_flux_x(ubar_R_plusx[i][j],ubar_L_plusx[i+1][j],gamma,c_h);
            }}

        //Update all real cells using x flux
        for(size_t i=2; i<num_xcells+2; i++){
            for(size_t j=2; j<num_ycells+2; j++){
                for(size_t k=0;k<9;k++){
                u[i][j][k]=u[i][j][k] - (dt/dx) * (flux_x[i][j][k]-flux_x[i-1][j][k]);
                }
            }
        }
  
        update_bcs(u,num_xcells,num_ycells);

        //Now perform y updates with the same timestep

        const auto& [ubar_L_plusy,ubar_R_plusy] = get_ubar_plus_yarrays(u,w,dt,dy,gamma,num_xcells,num_ycells,c_h);

    
         // update y flux array
        
        for(size_t i=1; i<num_xcells+2; i++){
            for(size_t j=1; j<num_ycells+2; j++){

                flux_y[i][j]=get_HLLC_flux_y(ubar_R_plusy[i][j],ubar_L_plusy[i][j+1],gamma,c_h);
            }}

        //Update all real cells using y flux
        for(size_t i=2; i<num_xcells+2; i++){
            for(size_t j=2; j<num_ycells+2; j++){
                for(size_t k=0;k<9;k++){
                u[i][j][k]=u[i][j][k] - (dt/dy) * (flux_y[i][j][k]-flux_y[i][j-1][k]);
                }
                
                //Now do the second partial source update for psi
               
                u[i][j][8]*=std::exp(-dt*(c_h*c_h)/c_p_squared);
         
            }
        }
        
        update_bcs(u,num_xcells,num_ycells);

      
        

   

    }while (t<tEnd);

}

    



