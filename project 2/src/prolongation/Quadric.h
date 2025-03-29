#ifndef _QUADRIC_
#define _QUADRIC_
#include"Prolongation.h"

template<int dim>
class Quadric:public Prolongation<dim>{
public:
    Vector operator()(const Vector &v) const;
};


#endif