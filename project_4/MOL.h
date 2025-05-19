#ifndef _MOL_
#define _MOL_
#include"Matrix.h"
#include"Function.h"
using namespace std;

class MOL{
protected:
    double r, t_begin, t_end, T, x_begin, x_end, k, h;
    int N, M; //M 表示空间离散的点数目， N表示时间的
    Matrix A;
    vector<Matrix> U;
public:
    MOL(const double &v, const double &h, const double &k, const double &t1, const double &t2, const double &x1, const double &x2);
    vector<vector<double>> get_bc(const bc_Function &);
    Matrix get_initial(const Cauchy_FUnction &);
    virtual void solve(const bc_Function &f, const Cauchy_FUnction &g)=0;
};

class FTCS:public MOL{
public:
    void solve(const bc_Function &f, const Cauchy_FUnction &g);
};


#endif