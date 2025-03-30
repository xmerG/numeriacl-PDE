#include"Multigrid.h"
#include"testFunction.h"
int main(){
    F1 f1;
    Multigrid<2> M1(f1, f1, BoundaryCondition::Dirichlet, 4);
    M1.print();
    return 0;
}