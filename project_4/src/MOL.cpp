#include"MOL.h"

MOL::MOL(const double &v, const double &_h, const double &_k, const double &t1, const double &t2, const double &x1, const double &x2):
            t_begin(t1), t_end(t2), x_begin(x1), x_end(x2), k(_k), h(_h){
    T=t2-t1;
    r=k*v/(h*h);
    M=(x_end-x_begin)/h-1;
    N=T/k;
    this->set_matrix();
}

void MOL::set_matrix(){
    for(int i=1; i<M; ++i){
        A.setValues(i, i, -2.0*r);
        A.setValues(i-1, i, r);
        A.setValues(i, i-1, r);
    }
    A.setValues(0,0, -2.0*r);
}

void MOL::get_bc(const bc_Function &bc){
    bc_value.resize(2, vector<double>(N+1, 0.0));
    for(int i=0; i<=N; ++i){
        bc_value[0][i]=bc(x_begin, k*i+t_begin);
        bc_value[1][i]=bc(x_end, k*i+t_begin);
    }
}

Vector MOL::get_initial(const Cauchy_FUnction & f){
    Vector initial(M);
    for(int i=0; i<M; ++i){
        initial.set_value(i,f(x_begin+h*(i+1)));
    }
    return move(initial);
}

void FTCS::solve(const bc_Function &f, const Cauchy_FUnction &g){
    Vector initial=move(MOL::get_initial(g));
    U.push_back(initial);
    this->get_bc(f);
    Triangle_Matrix I(M);
    for(int i=0; i<M; ++i){
        I.setValues(i, i, 1.0);
    }
    Triangle_Matrix B(M);
    B=move(I+A);
    for(int i=0; i<N; ++i){
        Vector b(M);
        b.set_value(0, bc_value[0][i]*r);
        b.set_value(M-1, bc_value[1][i]*r);
        Vector ui=move(B*initial+b);
        U.push_back(ui);
        initial=move(ui);
    }
}


void Crank_Nicolson::solve(const bc_Function &f, const Cauchy_FUnction &g){
    Vector initial=move(MOL::get_initial(g));
    U.push_back(initial);
    this->get_bc(f);
    Triangle_Matrix I(M);
    for(int i=0; i<M; ++i){
        I.setValues(i, i, 1.0);
    }
    Triangle_Matrix B1(M);
    Triangle_Matrix B2(M);
    B1=move(I+A*(-0.5));
    B2=move(I+A*0.5);
    for(int i=0; i<N; ++i){
        Vector b(M);
        b.set_value(0, (bc_value[0][i]+bc_value[0][i+1])*r*0.5);
        b.set_value(M-1, (bc_value[1][i]+bc_value[1][i+1])*0.5*r);
        Vector temp=move(B2*initial+b);
        Vector ui=move(B1/temp);
        U.push_back(ui);
        initial=move(ui);
    }
}

void BTCS::solve(const bc_Function &f, const Cauchy_FUnction &g){
    Vector initial=move(MOL::get_initial(g));
    U.push_back(initial);
    this->get_bc(f);
    Triangle_Matrix I(M);
    for(int i=0; i<M; ++i){
        I.setValues(i, i, 1.0);
    }
    Triangle_Matrix B(M);
    B=move(I-A);
    for(int i=0; i<N; ++i){
        Vector b(M);
        b.set_value(0, bc_value[0][i+1]*r);
        b.set_value(M-1, bc_value[1][i+1]*r);
        Vector temp=move(initial+b);
        Vector ui=move(B/temp);
        U.push_back(ui);
        initial=move(ui);
    }

}

void RegisterMOLMethods(){
    auto& factory = MOLFactory::getInstance();
    
    factory.registerMOL("FTCS", [](const double& v, const double& h, const double& k, 
                             const double& t1, const double& t2, 
                             const double& x1, const double& x2) {
    return std::make_unique<FTCS>(v, h, k, t1, t2, x1, x2);
    });
    
    factory.registerMOL("Crank_Nicolson", [](const double& v, const double& h, const double& k, 
                                           const double& t1, const double& t2, 
                                           const double& x1, const double& x2) {
        return std::make_unique<Crank_Nicolson>(v, h, k, t1, t2, x1, x2);
    });
    
    factory.registerMOL("BTCS", [](const double& v, const double& h, const double& k, 
                                 const double& t1, const double& t2, 
                                 const double& x1, const double& x2) {
        return std::make_unique<BTCS>(v, h, k, t1, t2, x1, x2);
    });
}