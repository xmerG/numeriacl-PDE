#include"Full_weighting.h"

template<int dim>
Vector Full_weighting<dim>::operator()(const Vector &v) const{
    if constexpr(dim==1){
        int n=v.getdim()+1;
        int currdim=n/2-1;
        Vector result(currdim);
        for(int i=0; i<currdim; ++i){
            double value=v(2*i)/4.0+v(2*i+1)/2.0+v(2*i+2)/4.0;
            result.set_Value(i, value);
        }
        return result;
    }
    else if constexpr(dim==2){
        
    }
}