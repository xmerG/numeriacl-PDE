#include"Quadratic.h"

template<int dim>
Vector Quadratic<dim>::operator()(const Vector &v) const{
    if constexpr(dim==1){
        int half_n=v.getdim()-1; //n/2
        int currdim=half_n*2+1;  //n+1
        Vector result(currdim);
        for(int i=0; i<=half_n; ++i){
            double value=(9*v(i)+9*v(i+1)-v(i-1)-v(i+2))/16.0;
            int newi=2*i;
            result.set_Value(newi, value);
            result.set_Value(newi+1, v(i));
        }
        double value= (3.0*v(0) + 6.0*v(1) -v(2))/8.0;
        result.set_Value(1, value);
        value=(3.0*v(half_n)+6.0*v(half_n-1)-v(half_n-2))/8.0;
        result.set_Value(currdim-2, value);
        return result;
    }
    else if constexpr(dim==2){
        int pre_dim=sqrt(v.getdim());   //n/2+1
        int half_n=pre_dim-1;  //n/2
        int resdim=2*half_n+1;   //n+1
        int quarter_n=half_n/2; //n/4
        Vector result(resdim*resdim);
        /*for(int i=0; i<half_n; ++i){
            for(int j=0; j<half_n; ++j){
                if(i==0 || i==half_n-1 || j==0 || j==half_n-1){
                    double value=0.25*(v(j, i)+v(j, i+1)+v(j+1, i)+v(j +1, i+1));
                    result.set_Value(2*j+1, 2*i+1,value);
                }
                else{
                    double value=0.5*(v(j, i+1)+v(j+1, i))+0.25*(v(j, i)+v(j+1, i+1))
                        -0.125*(v(j+1, i-1)+v(j+2, i)+v(j-1, i+1)+v(j, i+2));
                    result.set_Value(2*j+1, 2*i+1, value);
                }
            }
        }*/

        //拟合内部点
        /*for(int i=0; i<quarter_n; ++i){
            for(int j=0; j<quarter_n; ++j){
                double value=0.5*(v(j+1, i)+v(j, i+1))+0.25*v(j+1, i+1)-0.125*(v(j+2, i)+v(j, i+2));
                result.set_Value(2*j+1, 2*i+1, value);
                
                int indexj=j+quarter_n;
                value=0.5*(v(indexj, i)+v(indexj+1, i+1))+0.25*v(indexj, i+1)-0.125*(v(indexj-1, i)+v(indexj+1, i+2));
                result.set_Value(2*indexj+1, 2*i+1, value);

                int indexi=i+quarter_n;
                value=0.5*(v(j, indexi)+v(j+1, indexi+1))+0.25*v(j+1, indexi)-0.125*(v(j, indexi-1)+v(j+2, indexi+1));
                result.set_Value(2*j+1, 2*indexi+1, value);

                value=0.5*(v(indexj, indexi+1)+v(indexj+1, indexi))+0.25*v(indexj, indexi)
                        -0.125*(v(indexj-1, indexi+1)+v(indexj+1,indexi-1));
                result.set_Value(2*indexj+1, 2*indexi+1, value);
            }
        }*/
        for(int i=1; i<half_n-1; ++i){
            for(int j=1; j<half_n; ++j){
                double value=81.0*(v(j, i)+v(j+1, i+1)+v(j, i+1)+v(j+1, i))/256.0
                            -9.0*(v(j-1, i)+v(j-1, i+1)+v(j, i-1)+v(j+1, i-1)+v(j, i+2)+v(j+1, i+2)+v(j+2, i)+v(j+2, i+1))/256.0
                                    +(v(j-1, i-1)+v(j-1, i+2)+v(j+2, i+2)+v(j+2, i-1))/256.0;
                result.set_Value(2*j+1, 2*i+1, value);
            }
        }
        for(int i=1; i<half_n-1; ++i){
            double value=(-3.0*(v(0, i-1)+v(0, i+2))-6.0*(v(1, i-1)+v(1, i+2))+(v(2, i-1)+v(2, i+2))+
                            27.0*(v(0, i)+v(0, i+1))+54.0*(v(1, i)+v(1, i+1))/64.0-9.0*(v(2, i)+v(2, i+1)))/128.0;
            result.set_Value(1, 2*i+1, value);

            value=(-3.0*(v(half_n, i-1)+v(half_n, i+2))-6.0*(v(half_n-1, i-1)+v(half_n-1, i+2))+v(half_n-2, i-1)+v(half_n-2, i+2)+
                    27.0*(v(half_n, i)+v(half_n, i+1))+54.0*(v(half_n-1, i)+v(half_n-1, i+1))
                    -9.0*(v(half_n-2, i)+v(half_n-2, i+1)))/128.0;
            result.set_Value(resdim-2, 2*i+1, value);

            value=(-3.0*(v(i-1, 0)+v(i+2, 0))-6.0*(v(i-1, 1)+v(i+2, 1))+v(i-1, 2)+v(i+2, 2)+
                        27.0*(v(i, 0)+v(i+1, 0))+54.0*(v(i, 1)+v(i+1, 1))-9.0*(v(i, 2)+v(i+1, 2)))/128.0;
            result.set_Value(2*i+1, 1, value);
            
            value=(-3.0*(v(i-1, half_n)+v(i+2, half_n))-6.0*(v(i-1, half_n-1)+v(i+2, half_n-1))+v(i-1, half_n-2)+v(i+2, half_n-2)+
                        27.0*(v(i, half_n)+v(i+1, half_n))+54.0*(v(i, half_n-1)+v(i+1, half_n-1))
                        -9.0*(v(i, half_n-2)+v(i+1, half_n-2)))/128.0;
            result.set_Value(2*i+1, resdim-2, value);              
        }

        double value=(-3.0*(v(0, 2)+v(2,0))-6.0*(v(1, 2)+v(2, 1))+18.0*(v(0,1)+v(1,0))+9.0*v(0,0)+36.0*v(1,1)+v(2,2))/64.0;
        result.set_Value(1,1, value);

        value=(-3.0*(v(0, half_n-2)+v(2, half_n))-6.0*(v(1, half_n-2)+v(2, half_n-1))+18.0*(v(0, half_n-1)+v(1, half_n))
                +9.0*v(0, half_n)+36.0*v(1, half_n-1)+v(2, half_n-2))/64.0;
        result.set_Value(1, resdim-2);

        value=(-3.0*v(half_n-2, 0)+v(half_n, 2)-6.0*(v(half_n-2, 1)+v(half_n-1, 2))+18.0*(v(half_n-1, 0)+v(half_n, 1))
                +9.0*v(half_n, 0)+36.0*v(half_n-1, 1)+v(half_n-2, 2))/64.0;
        result.set_Value(resdim-2, 1, value);

        value=(-3.0*(v(half_n, half_n-2)+v(half_n-2, half_n))-6.0*(v(half_n-1, half_n-2)+v(half_n-2, half_n-1))+
                18.0*(v(half_n, half_n-1)+v(half_n-1, half_n))+9.0*v(half_n, half_n)+36.0*v(half_n-1, half_n-1)+
                v(half_n-2, half_n-2))/64.0;
        result.set_Value(resdim-2, resdim-2, value);

        //处理细网格中本来在粗网格上的点
        for(int i=0; i<=half_n; ++i){
            for(int j=0; j<=half_n; ++j){
                result.set_Value(2*j, 2*i, v(j, i));
            }
        }

        //处理边上的点
        for(int i=0; i<=half_n-1; ++i){
            for(int j=0; j<=half_n; ++j){
                double value=(9.0*v(j, i)+9.0*v(j, i+1)-v(j, i-1)-v(j, i+2))/16.0;
                result.set_Value(2*j, 2*i+1, value);

                value=(9.0*v(j, i)+9.0*v(j+1, i)-v(j-1, i)-v(j+2, i))/16.0;
                result.set_Value(2*j+1, 2*i, value);
            }
        }


        for(int i=0; i<=half_n; ++i){
            double value=(3.0*v(i, 0)+6.0*v(i, 1)-v(i, 2))/8.0;
            result.set_Value(2*i, 1, value);

            value=(3.0*v(i, half_n)+6.0*v(i, half_n-1)-v(i, half_n-2))/8.0;
            result.set_Value(2*i, resdim-2, value);

            value=(3.0*v(0, i)+6.0*v(1, i)-v(2, i))/8.0;
            result.set_Value(1, 2*i, value);

            value=(3.0*v(half_n, i)+6.0*v(half_n-1, i)-v(half_n-2, i))/8.0;
            result.set_Value(resdim-2, 2*i, value);
            
        }
        /*for(int i=0; i<=half_n; ++i){
            for(int j=0; j<=half_n; ++j){
                int newi=2*i;
                int newj=2*j;
                result.set_Value(newj, newi, v(j,i));

                double value=(9*v(j, i)+9*v(j+1, i)-v(j-1, i)- v(j+2, i))/16.0;
                result.set_Value(newj+1, newi, value);

                value=(9*v(j, i)+9*v(j, i+1)-v(j, i-1)-v(j,i+2))/16.0;
                result.set_Value(newj, newi+1, value);

                value=0.5*(v(j, i+1)+v(j+1, i))+0.25*v(j+1,i+1)-0.125*(v(j, i+2)+v(j+2, i));
                result.set_Value(newj, newi, value);
            }
            //---------------------------------to be modified----------------------------------------------------------------
        }


        int index=half_n-1;
        double value=0.5*(v(index, half_n)+v(half_n, index))+0.25*(v(index, index))-0.125*(v(index-1, half_n)+v(half_n, index-1));
        result.set_Value(resdim-2, resdim-2, value); */
        return result;
    }
    else{
        cerr<<"can't solve when dimension is over 2!"<<endl;
        return Vector();
    }
}

template class Quadratic<1>;
template class Quadratic<2>;