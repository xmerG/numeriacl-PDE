#ifndef _RESTRICTION_
#define _RESTRICTION_
#include"Vector.h"

template<int dim>
class Retriction{
public:
    virtual Vector operator()(const Vector &v) const=0;
};


#endif