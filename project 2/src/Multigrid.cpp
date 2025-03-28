#include"Multigrid.h"

template<int dim>
Vector Multigrid<dim>::w_Jacobi(const int &i){
    Sparse_Matrix A=discretors[i].first;
    Vector v=discretors[i].second;
    const double weighted=1.0/3.0;
    double tempcoeff=weighted*pow((i*h),2);
    Vector temp=A*tempcoeff;
    return v-temp*v;
}
