#include"Multigrid.h"
#include"testFunction.h"
int main(){
    F1 f1;
    /*Neumann g1;
    Laplacian f;
    int n=256;
    Vector v(n-1);
    Vector &v0=v;
    double value=f1(1.0/n);
    Multigrid<1> M1(f, g1, BoundaryCondition::Neumann, n, value);
    M1.solve("injection", "linear","v-cylce", v0, 5,5);
    M1.print_to_file("output.json", BoundaryCondition::Neumann, f1);

    Multigrid<1> M3(f, f1, BoundaryCondition::Dirichlet, n);
    M3.solve("injection", "linear","v-cylce", v0, 5,5);
    M3.print_to_file("output.json", BoundaryCondition::Dirichlet, f1);*/

    int n=8;
    F2 primitive;
    Laplacian2 l;
    Vector V((n-1)*(n-1));
    Vector &v1=V;
    Multigrid<2> M2(l, primitive, BoundaryCondition::Dirichlet, n);
    M2.solve("injection", "linear","v-cylce", v1, 5,5);
    M2.print_to_file("output.json", BoundaryCondition::Dirichlet, primitive);



    return 0;
}