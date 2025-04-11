#ifndef _QUADRATIC_
#define _QUADRATIC_
#include"Prolongation.h"

template<int dim>
class Quadratic:public Prolongation<dim>{
public:
    Vector operator()(const Vector &v) const;
};


#endif