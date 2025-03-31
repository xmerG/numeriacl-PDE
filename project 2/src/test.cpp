#include"Multigrid.h"
#include"testFunction.h"
int main(){
    F1 f1;
    Laplacian f;
    Vector v(7);
    Vector &v0=v;
    Multigrid<1> M1(f, f1, BoundaryCondition::Dirichlet, 8);
    M1.solve("full_weighting", "linear","FMG", v0, 2,2);
    M1.print_to_file("output.json", BoundaryCondition::Dirichlet, f1);
    return 0;
}