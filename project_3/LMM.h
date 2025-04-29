#ifndef _LMM_
#define _LMM_

#include"IVP_solver.h"

class LMM:public IVP_solver{
protected:
    int step; //step
    Matrix alpha=Matrix(1, step);
    Matrix beta=Matrix(1, step);
public:
    virtual void OneStep(const Function &f, int n)=0;
    void setInitials(const Function &f);
    void solve(const Function &f);
    void solve(const Function &f, const double &_initial, const double &, const double &, int);
};


class Adams_Bashforth:public LMM{
public:
    Adams_Bashforth();
    Adams_Bashforth(int p);
    void OneStep(const Function &f, const Matrix &);
};

#endif