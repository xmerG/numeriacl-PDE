#ifndef _FULL_WEIGHTING_
#define _FULL_WEIGHTING_
#include"Restriction.h"

template<int dim>
class Full_weighting:public Retriction<dim>{
public:
    Vector operator()(const Vector &v) const;
};



#endif