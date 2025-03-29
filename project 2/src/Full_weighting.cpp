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
        int n=sqrt(v.getdim())+1;
        int temp=n/2-1;
        int currdim=pow(temp,2);
        Vector result(currdim);
        for(int i=0; i<temp; ++i){
            for(int j=0; j<temp; ++j){
                int newi=2*i*(n-1);
                int newj=2*j;
                double value=v(newi+newj)+2*v(newi+newj+1)+v(newi+newj+2);
                newi+=n-1;
                value+=2*v(newi+newj)+4*v(newi+newj+1)+2*v(newi+newj+2);
                newi+=n-1;
                value+=v(newi+newj)+2*v(newi+newj+1)+v(newi+newj+2);
                value=value/16.0;
                result.set_Value(temp*i+j, value);
            }
        }
        return result;
    }
    else{
        cerr<<"can't solve problems over 2 dimension!"<<endl;
        return Vector();
    }
}

template class Full_weighting<1>;
template class Full_weighting<2>;