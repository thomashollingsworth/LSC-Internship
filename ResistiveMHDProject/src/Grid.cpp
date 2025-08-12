#include "Grid.h"


Grid::Grid(size_t nx, size_t ny,size_t g_)
    :ghost_cells(g_),
    num_xcells(nx),
    num_ycells(ny),
    

    x(nx),
    y(ny),

    U(nx,ny,g_),
    Prim(nx,ny,g_),
    uBarL(nx,ny,g_),
    uBarR(nx,ny,g_),
    primL(nx,ny,g_),
    primR(nx,ny,g_),
    Delta(nx,ny,g_),
    Chi(nx,ny,g_),
    fluxL(nx,ny,g_),
    fluxR(nx,ny,g_),
    fluxHLLD(nx,ny,g_),
    c_fx(nx,ny,g_),
    c_fy(nx,ny,g_),
    c_fz(nx,ny,g_),
    uStarL(nx,ny,g_),
    uStarR(nx,ny,g_),
    uDoubleStar(nx,ny,g_)

     {};

void Grid::reset_intermediate_arrays(){
 
    for(size_t i=0;i<num_xcells+2*ghost_cells;i++){
        for(size_t j=0;j<num_ycells+2*ghost_cells;j++){
        uBarL(i,j)*=0;
        uBarR(i,j)*=0;
        primL(i,j)*=0;
        primR(i,j)*=0;
        Delta(i,j)*=0;
        Chi(i,j)*=0;
        fluxL(i,j)*=0;
        fluxR(i,j)*=0;
        fluxHLLD(i,j)*=0;
        c_fx(i,j)*=0;
        c_fy(i,j)*=0;
        c_fz(i,j)*=0;
        uStarL(i,j)*=0;
        uStarR(i,j)*=0;
        uDoubleStar(i,j)*=0;

}}
}

