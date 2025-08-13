#include "SavingRoutine.h"


void save_to_file(const Grid& grid, const SimulationConfig& cfg, size_t iteration,double time){
    if(cfg.save){
        if(iteration%cfg.save_interval==0){
            const size_t nx = grid.num_xcells;
            const size_t ny = grid.num_ycells;
            const size_t g  = grid.ghost_cells;

             //Output the data
            std::filesystem::create_directories(cfg.save_directory);

            std::string filename= "t_" + std::to_string(time);
            filename=cfg.save_directory+ "/" + filename;

            std::ofstream outfile(filename);

            for(size_t i=0; i<nx;i++){
                for(size_t j=0; j<ny;j++){
            outfile << grid.x[i] << " " << grid.y[j] << " " << grid.Prim(i+g,j+g).density() << " "<< grid.Prim(i+g,j+g).velocity().x()<< " "<< grid.Prim(i+g,j+g).velocity().y()<< " "<< grid.Prim(i+g,j+g).velocity().z()<< " "<< grid.Prim(i+g,j+g).pressure()<< " "<< grid.Prim(i+g,j+g).B().x()<< " "<< grid.Prim(i+g,j+g).B().y()<< " "<< grid.Prim(i+g,j+g).B().z()<< " "<< grid.Prim(i+g,j+g).psi()<<"\n";
            }
            outfile<<"\n"; //Blank line between y rows for correct gnuplot formatting
                }
        }
    }

    if(cfg.plot){
        if(iteration%cfg.plot_interval==0){
             //Output the data
            std::string filename= "t_" + std::to_string(time);
            std::filesystem::create_directories(cfg.save_directory);
            filename=cfg.save_directory+"/"+filename;
            const size_t nx = grid.num_xcells;
            const size_t ny = grid.num_ycells;
            const size_t g  = grid.ghost_cells;
             
             if(cfg.save_interval%iteration!=0){ 
                std::ofstream outfile(filename);

                for(size_t i=0; i<nx;i++){
                for(size_t j=0; j<ny;j++){
                outfile << grid.x[i] << " " << grid.y[j] << " " << grid.Prim(i+g,j+g).density() << " "<< grid.Prim(i+g,j+g).velocity().x()<< " "<< grid.Prim(i+g,j+g).velocity().y()<< " "<< grid.Prim(i+g,j+g).velocity().z()<< " "<< grid.Prim(i+g,j+g).pressure()<< " "<< grid.Prim(i+g,j+g).B().x()<< " "<< grid.Prim(i+g,j+g).B().y()<< " "<< grid.Prim(i+g,j+g).B().z()<< " "<< grid.Prim(i+g,j+g).psi()<<"\n";
                }
                outfile<<"\n"; //Blank line between y rows for correct gnuplot formatting
                }}
            
            std::string gnuplot_cmd = "gnuplot -e \"filename='" + filename + "'\" plot/plot_script.plt";
            std::system(gnuplot_cmd.c_str());
            std::cout << "Plot saved to " <<filename<< ".png" << std::endl;
        }}
    }