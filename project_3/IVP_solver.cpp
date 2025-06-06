#include"IVP_solver.h"

IVP_solver::IVP_solver(){}

IVP_solver::IVP_solver(Method _m, const int &_d):method(_m), d(_d){}

IVP_solver::IVP_solver(Method m, const int & _d, const Matrix &_initial, const double & t1, const double & t2, int N):
            method(m), d(_d), t_begin(t1), t_end(t2), maxStep(N){
                k=(t_end-t_begin)/N;
                solution.push_back(_initial);
            }

void IVP_solver::setValues(const double & t1, const double & t2, int N){
    t_begin=t1;
    t_end=t2;
    maxStep=N;
    k=(t_end-t_begin)/N;
}


Matrix IVP_solver::getSolution(const int &index){
    Matrix U_new;
    U_new=solution[index];
    return move(U_new);
}
