#ifndef _INJECTION_
#define _INJECTION_
#include"Restriction.h"

template<int dim>
class Injection:public Retriction<dim>{
public:
    Vector operator()(const Vector &v) const;
};

#endif