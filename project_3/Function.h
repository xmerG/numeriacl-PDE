#ifndef _FUNCTION_
#define _FUNCTION_
#include"Matrix.h"

class Function{
public:
    virtual Matrix operator()(const Matrix &u, const double &t) const=0;
};


#endif