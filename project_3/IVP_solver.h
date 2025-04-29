#ifndef _IVP_SOLVER_
#define _IVP_SOLVER_
#include"Matrix.h"
#include"Method.h"
#include"Function.h"

class IVP_solver{
protected:
    int d; //dimension of the ode
    double t_begin, t_end; //solve the equation on [t_begin, t_end]
    double k; //size of one step
    int current; //current n
    bool explicit_method=true;
    //Matrix initial=move(Matrix(d, 1));
    int maxStep; 
    int accuracy=1; // accuracy 
    Matrix solution=move(Matrix(1, d));
    Method method;
public:
    IVP_solver();
    IVP_solver(Method, const int &);
    IVP_solver(Method, const int &, const Matrix &_initial, const double &, const double &, int);
    void setValues(const double &, const double &, int);
    virtual void solve(const Function &f)=0;
    virtual void solve(const Function &f, const Matrix &_initial, const double &, const double &, int)=0;
};

#endif