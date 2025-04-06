#include"Full_weighting.h"

template<int dim>
Vector Full_weighting<dim>::operator()(const Vector &v) const{
    if constexpr(dim==1){
        int n=v.getdim()-1;  //n
        int currdim=n/2+1;  //n/2+1
        int half_n=currdim-1;
        Vector result(currdim);
        for(int i=1; i<half_n; ++i){
            double value=v(2*i-1)/4.0+v(2*i)/2.0+v(2*i+1)/4.0;
            result.set_Value(i, value);
        }
        result.set_Value(0, v(0));
        result.set_Value(half_n, v(n));
        return result;
    }
    else if constexpr(dim==2){
        int n=sqrt(v.getdim())-1;  //n
        int temp=n/2+1;  //n/2+1
        int currdim=pow(temp,2);  //(n/2+1)*(n/2+1)
        int half_n=temp-1;
        Vector result(currdim);
        for(int i=1; i<half_n; ++i){
            for(int j=1; j<half_n; ++j){
                int newi=2*i;
                int newj=2*j;
                double value=4*v(newj,newi)+2*(v(newj,newi-1)+v(newj, newi+1)+v(newj-1,newi)+v(newj+1, newi))
                        +v(newj-1, newi-1)+v(newj-1, newi+1)+v(newj+1, newi-1)+v(newj+1, newi+1);
                value/=16.0;
                result.set_Value(j,i, value);
            }
        }
        for(int i=0; i<=half_n; ++i){
            int newi=2*i;
            result.set_Value(0, i, v(0, newi));
            result.set_Value(half_n, v(n,  newi));
            result.set_Value(i, 0, v(newi, 0));
            result.set_Value(i, half_n, v(newi, n));
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