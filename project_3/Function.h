#ifndef _FUNCTION_
#define _FUNCTION_
#include"Matrix.h"

class Function{
public:
    virtual Matrix operator()(const Matrix &u, const double &t) const=0;
    Matrix df_x(const Matrix &u, const double &t){
        int dim=u.get_col();
        Matrix D(dim+1, dim+1);
        double delta=1e-7;
        Matrix u1;
        u1=u;
        Matrix u2;
        u2=u;
        for(int i=0; i<=dim; ++i){
            u1.set_elements(1, i, u1(1, i)+delta);
            u2.set_elements(1, i, u2(1, i)-delta);
            Matrix f_U_i=move((this->operator()(u1, t)-this->operator()(u2, t))*(1.0/delta));
            D.set_elements(f_U_i, i, i, 0, dim-1);
        }
        return move(D);
    }
};


#endif