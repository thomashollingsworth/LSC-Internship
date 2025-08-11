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
    Delta(nx,ny,g_),
    Chi(nx,ny,g_),
    fluxL(nx,ny,g_),
    fluxR(nx,ny,g_),
    c_fx(nx,ny,g_),
    c_fy(nx,ny,g_),
    c_fz(nx,ny,g_),
    S_L(nx,ny,g_),
    S_R(nx,ny,g_),
    S_M(nx,ny,g_)

     {};


