#ifndef _LINEAR_
#define _LINEAR_
#include"Prolongation.h"

template<int dim>
class Linear:public Prolongation<dim>{
public:
    Vector operator()(const Vector &v) const;
};


#endif