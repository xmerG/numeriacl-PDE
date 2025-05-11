#ifndef _RUNGE_KUTTA_
#define _RUNGE_KUTTA_
#include"IVP_solver.h"

class Implicit_Runge_Kutta:public IVP_solver{
protected:
    int s;
    vector<vector<double>> a=vector<vector<double>>(s, vector<double>(s, 0.0));
    vector<double> b=vector<double>(s, 0.0);
    vector<double> c=vector<double>(s, 0.0);
    void OneStep(const Function &f, int n);
    void solve(const Function &f);
};

class classical_RK:public IVP_solver{
private:
    void OneStep(const Function &f, int n);
};

class ESDIRK:public Implicit_Runge_Kutta{

};




#endif