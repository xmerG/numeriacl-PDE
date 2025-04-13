#include"Multigrid.h"
#include"testFunction.h"
#include <chrono>
int main(){
    vector<double> corsa{0.0, 1.0, 3.0, 5.0, 7.0, 4.0, 8.0, 6.0, 10.0};
    Quadratic<2> q;
    Vector corsa_v(corsa);
    Vector fine=q(corsa_v);
    fine.print();

    F1 f1;
    //Neumann g1;
    //Laplacian f;
    int n=128;
    Vector v(n+1);
    Vector &v0=v;
    double value=f1(0.0);
    //Multigrid<2> M1(f, g1, BoundaryCondition::Mixed, n, vector<int>{0,1,0,1});
    cout<<"------------------------------------------------------------------"<<endl;
    //Multigrid<1> M1(f, g1, BoundaryCondition::Neumann, n, vector<int>{0,1});
    //M1.solve("full_weighting", "quadric","FMG", v0, 5,5, 1e-8, value);
    //M1.print_to_file("output.json", f1);

    //Multigrid<1> M3(f, f1, BoundaryCondition::Dirichlet, n);
    //M3.solve("injection", "linear","FMG", v0, 5,5,1e-8);
    //M3.print_to_file("output.json", f1);

    //int n=8;
    F2 primitive;
    Laplacian2 l;
    Vector V((n+1)*(n+1));
    Vector &v1=V;
    //Multigrid<2> M2(l, primitive, BoundaryCondition::Dirichlet, n);
    //M2.solve("full_weighting", "quadric","FMG", v1, 5,5, 10e-8);
    //M2.print_to_file("output.json", primitive);

    Neumann2 g2;
    value=primitive(0.0, 0.0);
    Multigrid<2> M3(l, g2, BoundaryCondition::Neumann, n);
    M3.solve("full_weighting", "quadric","FMG", v1, 5,5, 1e-8, value);
    M3.print_to_file("output.json", primitive);

    vector<int> mixed=vector<int>{0,1,1,0};
    Mixed2 m(mixed);
    Multigrid<2> M4(l, m, BoundaryCondition::Mixed, n, mixed);
    M4.solve("full_weighting", "linear","FMG", v1, 5,5, 10e-8);
    M4.print_to_file("output.json", primitive);



    return 0;
}