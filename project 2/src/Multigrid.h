#ifndef _MULTIGRID_
#define _MULTIGRID_
#include"Sparse_Matrix.h"
#include<cmath>

using couplet=pair<Sparse_Matrix&, Vector&>;

template<int dim>
class Multigrid{
private:
    map<int, couplet> discretors; //int 代表网格大小，如2代表2h的网格 
    double h; //size of a single grid
    int n; //number of grids
    Vector w_Jacobi(const int &i);
    Vector (*Restriction)(int i, Vector &v);
};


#endif