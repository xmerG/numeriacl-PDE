#ifndef _MULTIGRID_
#define _MULTIGRID_
#include"Sparse_Matrix.h"
#include"restriction/Full_weighting.h"
#include"restriction/Injection.h"
#include"prolongation/Linear.h"
#include"prolongation/Quadric.h"
#include"BoudaryCondition.h"
#include"Function.hpp"
#include <memory>
#include <lapacke.h>


using couplet=pair<Sparse_Matrix, Vector>;

template<int dim>
class Multigrid{
private:
    map<int, couplet> discretors; //int 代表网格数目，如2代表2的网格数，4代表4的网格数 
    int n;
    Vector solutions;
    unique_ptr<Retriction<dim>> restriction;
    unique_ptr<Prolongation<dim>> prolongation;
    Vector w_Jacobi(int i, const Vector &initial);
    Vector(Multigrid<dim>::*cycle)(const int &n, Vector& initial_guess, int nu1, int nu2);
    Vector V_cycle(const int &n, Vector& initial_guess, int nu1, int nu2);
    Vector FMG(const int &n, Vector& initial_guess, int nu1, int nu2);
    Vector corsa_solve(const int &i);
    void create_grids_D(const Function &f, const Function &g,const int &i);
    void create_grids_N(const Function &f, const Function &g,const int &i, double value);
    vector<double> error(const Function &f);

public:
    Multigrid();
    Multigrid(const Function &f, const Function &g, BoundaryCondition bc, const int &i);
    void solve(string restriction, string prolongation, string cycle, Vector& initial_guess, int nu1, int nu2);
    void print();
    void print_to_file(const string &filename, BoundaryCondition BC,const Function &f);
};


#endif