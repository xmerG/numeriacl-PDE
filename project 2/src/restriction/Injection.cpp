#include"Injection.h"
template<int dim>
Vector Injection<dim>::operator()(const Vector &v) const{
    if constexpr (dim == 1){
        int n=v.getdim()-1;  //n
        int currdim=n/2+1;   //n/2+1
        Vector result(currdim);
        for(int j=0; j<currdim; ++j){
            result.set_Value(j, v(2*j));
        }
        return result;
    }
    else if constexpr(dim==2){
        int n=sqrt(v.getdim())-1; //n/2
        int temp=n/2+1;  //n/2+1
        int currdim=pow(temp,2);
        Vector result(currdim);
        for(int i=0; i<temp; ++i){
            for(int j=0; j<temp; ++j){
                result.set_Value(j, i, v(2*j, 2*i));
            }
        }
        return result;
    }
    else{
        std::cerr<<"can't solve problems over dimension 3 !"<<std::endl;
        return Vector();
    }
}
template class Injection<1>;
template class Injection<2>;