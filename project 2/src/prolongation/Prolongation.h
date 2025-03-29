#ifndef _PROLONGATION_
#define _PROLONGATION_

#include"../Vector.h"

template<int dim>
class Prolongation{
public:
    virtual Vector operator()(const Vector &v) const=0;
};


#endif