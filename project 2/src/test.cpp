#include"prolongation/Linear.h"

int main(){
    vector<double> e{1.0, 7.0, 5.0, 4.0, 2.0, 0.0, -1.0, 3.0, 8.0};
    Vector v1=Vector(e);
    Linear<2> inj1;
    Vector v2=inj1(v1);
    v2.print();
    return 0;
}