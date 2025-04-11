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



using couplet=pair<Sparse_Matrix, Vector>;

template<int dim>
class Multigrid{
private:
    map<int, couplet> discretors; //int 代表网格数目，如2代表2的网格数，4代表4的网格数 
    int n;  //网格个数
    Vector solutions;
    BoundaryCondition BC;
    int counter=1; //迭代次数 
    unique_ptr<Retriction<dim>> restriction;
    unique_ptr<Prolongation<dim>> prolongation;
    Vector w_Jacobi(int i, const Vector &initial);
    Vector V_cycle(const int &n, Vector& initial_guess, int nu1, int nu2);
    Vector FMG(const int &_n, int nu1, int nu2);
    vector<double> corsa_solve(const int &i);
    void create_grids_D(const Function &f, const Function &g,const int &i);
    void create_grids_N(const Function &f, const Function &g,const int &i);
    void create_grids_M(const Function &f, const Function &g,const int &i, 
        const vector<int> &mixed=vector<int> {0,0,0,0});
    Vector error(const Function &f);

public:
    Multigrid();
    Multigrid(const Function &f, const Function &g, BoundaryCondition bc, const int &i,
            const vector<int> &mixed=vector<int>{0,0,0,0});
    void solve(const string &restriction, const string &prolongation, const string &cycle, Vector& initial_guess, int nu1, int nu2,
        double tol,const double &value=0.0, int max_itr=50);
    void print();
    void print_to_file(const string &filename, const Function &f);
};


#endif