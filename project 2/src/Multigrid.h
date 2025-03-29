#ifndef _MULTIGRID_
#define _MULTIGRID_
#include"Sparse_Matrix.h"
#include"restriction/Full_weighting.h"
#include"restriction/Injection.h"
#include"prolongation/Linear.h"
#include"prolongation/Quadric.h"
#include <memory>

using couplet=pair<Sparse_Matrix&, Vector&>;

template<int dim>
class Multigrid{
private:
    map<int, couplet> discretors; //int 代表网格大小，如2代表2h的网格 
    double h; //size of a single grid
    int n; //number of grids
    unique_ptr<Retriction<dim>> restriction;
    unique_ptr<Prolongation<dim>> prolongation;
    Vector w_Jacobi(const int &i);
    Vector 2_grid_correction(Vector &initial_guess, Vector &f, int size, int v1, int v2);
    Vector V_cycle();
    Vector FMG();
};


#endif