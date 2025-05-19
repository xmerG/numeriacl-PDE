#include"MOL.h"

MOL::MOL(const double &v, const double &_h, const double &_k, const double &t1, const double &t2, const double &x1, const double &x2):
            t_begin(t1), t_end(t2), x_begin(x1), x_end(x2), k(_k), h(_h){
    T=t2-t1;
    r=k*v/(h*h);
    M=(x_end-x_begin)/h-1;
    N=T/k;
}

vector<vector<double>> MOL::get_bc(const bc_Function &bc){
    vector<vector<double>> boundary_value(2, vector<double>(N+1, 0.0));
    for(int i=0; i<=N; ++i){
        boundary_value[0][i]=bc(x_begin, k*i+t_begin);
        boundary_value[1][i]=bc(x_end, k*i+t_begin);
    }
    return boundary_value;
}

Matrix MOL::get_initial(const Cauchy_FUnction & f){
    Matrix initial(M,1);
    for(int i=0; i<M; ++i){
        initial.set_elements(i, 0, f(x_begin+h*(i+1)));
    }
    return move(initial);
}

void FTCS::solve(const bc_Function &f, const Cauchy_FUnction &g){
    Matrix initial=move(MOL::get_initial(g));
    Matrix b(M, 1);
    
}