#include"Quadric.h"

template<int dim>
Vector Quadric<dim>::operator()(const Vector &v) const{
    if constexpr(dim==1){
        int half_n=v.getdim()+1; //n/2
        int currdim=half_n*2-1;  //n-1
        Vector result(currdim);
        for(int i=0; i<half_n; ++i){
            double value=(9*v(i)+9*v(i-1)-v(i-2)-v(i+1))/16.0;
            int newi=2*i;
            result.set_Value(newi, value);
            result.set_Value(newi+1, v(i));
        }
        return result;

    }
    else if constexpr(dim==2){
        int pre_dim=sqrt(v.getdim());   //n/2-1
        int half_n=pre_dim+1;  //n/2
        int resdim=2*half_n-1;   //n-1
        Vector result(resdim*resdim);
        for(int i=0; i<half_n; ++i){
            for(int j=0; j<half_n; ++j){
                int newi=2*i;
                int newj=2*j;
                result.set_Value(newi+1, newj+1, v(j,i));

                double value=(9*v(j, i)+9*v(j, i-1)-v(j, i-2)- v(j, i+1))/16.0;
                result.set_Value(newj+1, newi, value);

                value=(9*v(j, i)+9*v(j-1, i)-v(j+1, i)-v(j-2,i))/16.0;
                result.set_Value(newj, newi+1, value);

                value=0.5*(v(j, i-1)+v(j-1, i))+0.25*v(j,i)-0.125*(v(j+1, i-1)+v(j-1, i+1));
                result.set_Value(newj, newi, value);
            }
        }
        return result;
    }
    else{
        cerr<<"can't solve when dimension is over 2!"<<endl;
        return Vector();
    }
}

template class Quadric<1>;
template class Quadric<2>;