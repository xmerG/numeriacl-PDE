#ifndef _LMM_
#define _LMM_

#include"IVP_solver.h"

class LMM:public IVP_solver{
protected:
    int step; //step
    vector<double> alpha=vector<double>(step, 0.0);  //(a1, a2,..., as-1)
    vector<double> beta=vector<double>(step+1, 0.0);
    vector<Matrix> initials;
    virtual void OneStep(const Function &f, int n)=0;
public:
    void setInitials(const Function &f);
    void solve(const Function &f);
    void solve(const Function &f, const double &_initial, const double &, const double &, int);
    void solve(const Function &f, const double &_initial, const double &, const double &, int, const int &);
};


class Adams_Bashforth:public LMM{
private:
    void set_beta(int p);
    void OneStep(const Function &f,  int n);
public:
    Adams_Bashforth();
    Adams_Bashforth(int p);

};

class Adams_Moulton:public LMM{
private:
    void set_beta(int p);
    void OneStep(const Function &, int);
public:
    Adams_Moulton();
    Adams_Moulton(int p);
};

class Backward_Differential:public LMM{
private:
    void set(int p);
    void OneStep(const Function &f, int n);
public:
    Backward_Differential();
    Backward_Differential(int p);
};

#endif