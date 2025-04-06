#include"Linear.h"

template<int dim>
Vector Linear<dim>::operator()(const Vector &v) const{
    if constexpr(dim==1){
        int half_n=v.getdim()-1;  //n/2
        int resdim=half_n*2+1;  //n+1
        Vector result(resdim);
        for(int j=0; j<=half_n; ++j){
            double value=0.5*v(j)+0.5*v(j+1);
            result.set_Value(2*j, v(j));
            result.set_Value(2*j+1, value);
        }
        return result;
    }
    else if constexpr(dim==2){
        int pre_dim=sqrt(v.getdim());   //n/2+1
        int half_n=pre_dim-1;  //n/2
        int resdim=2*half_n+1;   //n+1
        Vector result(resdim*resdim);
        for(int i=0; i<=half_n; ++i){
            for(int j=0; j<=half_n; ++j){ 
                int newi=2*i;
                int newj=2*j;               
                result.set_Value(newj,newi, v(j, i));
                
                double value=0.5*(v(j,i)+v(j+1, i));
                result.set_Value(newj+1,newi, value);

                value=0.5*(v(j,i)+v(j,i+1));
                result.set_Value(newj,newi+1, value);

                value=value*0.5+0.25*(v(j+1, i)+v(j+1, i+1));
                result.set_Value(newj+1, newi+1,value);
            }
        }
        return result;
    }
    else{
        cerr<<"can't solve when dimension is over 2 !"<<endl;
        return Vector();
    }
}

template class Linear<1>;
template class Linear<2>;