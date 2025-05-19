#ifndef _FUNCTION_
#define _FUNCTION_
#include<cmath>

class bc_Function{
public:
    virtual double operator()(const double &x, const double &y)const=0;
};

class Cauchy_FUnction{
public:
    virtual double operator()(const double &x) const=0;
};



#endif