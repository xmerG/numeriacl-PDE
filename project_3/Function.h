#ifndef _FUNCTION_
#define _FUNCTION_
#include"Matrix.h"

class Function{
public:
    virtual double operator()(const Matrix &u, const double &t)=0;
};


#endif