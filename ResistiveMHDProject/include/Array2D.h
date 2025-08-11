#ifndef ARRAY2D_H
#define ARRAY2D_H
#include <vector>

//This structure will be a general template for a std::vector that represents a 2D array

template<typename Component>
struct Array2D{
   
    size_t nx;
    size_t ny;
    size_t ghost_cells;
    std::vector<Component> data;
   
    //Constructor
    Array2D(size_t nx,size_t ny,size_t ghost_cells=2)
    : nx(nx), ny(ny), ghost_cells(ghost_cells), data((nx+2*ghost_cells)*(ny+2*ghost_cells)) {}
    
    
    //Access indices by 2D indexing
    Component& operator()(size_t i, size_t j){
    return data[(ny+2*ghost_cells)*i+j];}
    
    const Component& operator()(size_t i, size_t j) const {
    return data[(ny+2*ghost_cells)*i+j];}

};



#endif