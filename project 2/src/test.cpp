#include"Multigrid.h"
#include"testFunction.h"
#include <chrono>
int main(){
    //vector<double> fine{0, 2.5, 5, 4.5, 4, 1, 4, 7, 6.5, 6, 2, 5.5, 9, 8.5, 8, 4.5, 4.75, 5, 5.25, 5.5, 7, 4, 1, 2, 3 };
    //Injection<2> pro;
    //Vector fine_v(fine);
    //Vector corsa=pro(fine_v);
    //corsa.print();

    F1 f1;
    Neumann g1;
    Laplacian f;
    int n=4;
    Vector v(n+1);
    Vector &v0=v;
    double value=f1(0.0);
    Multigrid<2> M1(f, g1, BoundaryCondition::Mixed, n, vector<int>{0,1,0,1});
    cout<<"------------------------------------------------------------------"<<endl;
    Multigrid<1> M2(f, g1, BoundaryCondition::Neumann, n, vector<int>{0,1});
    //M1.solve("injection", "linear","FMG", v0, 5,5, 1e-8, value);
    //M1.print_to_file("output.json", f1);

    //Multigrid<1> M3(f, f1, BoundaryCondition::Dirichlet, n);
    //M3.solve("injection", "linear","FMG", v0, 5,5,1e-8);
    //M3.print_to_file("output.json", f1);

    //int n=8;
    //F2 primitive;
    //Laplacian2 l;
    //Vector V((n+1)*(n+1));
    //Vector &v1=V;
    //Multigrid<2> M2(l, primitive, BoundaryCondition::Dirichlet, n);
    //M2.solve("injection", "linear","v-cylce", v1, 5,5, 10e-8);
    //M2.print_to_file("output.json", primitive);

    //Neumann2 g2;
    //value=primitive(0.0, 0.0);
    //Multigrid<2> M3(l, g2, BoundaryCondition::Neumann, n);
    //M3.solve("injection", "linear","FMG", v1, 5,5, 10e-8, value);
    //M3.print_to_file("output.json", primitive);



    return 0;
}