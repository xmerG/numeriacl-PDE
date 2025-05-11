#include"LMM.h"

void LMM::solve(const Function &f){
    this->setInitials(f);
    for(int i=0; i<=maxStep-step; ++i){
        this->OneStep(f, i);
    }
}

void LMM::solve(const Function &f, const double &_initial, const double &t1, const double &t2, int N){
    IVP_solver::setValues(t1, t2, N);
    LMM::solve(f);
}


void LMM::setInitials(const Function &f){
    if(step>1){
        
    }
}


Adams_Bashforth::Adams_Bashforth(){
    method=Method::ABF;
    this->set_beta(1);
}

Adams_Bashforth::Adams_Bashforth(int p){
    method=Method::ABF;
    accuracy=p;
    this->set_beta(accuracy);
}

void Adams_Bashforth::OneStep(const Function &f, int n){
    Matrix U_new(1, d);
    for(int j=0; j<step; ++j){
        U_new=move(U_new+f(solution[n+j], k*(n+j))*(beta[j]*k+t_begin));
    }
    U_new=move(U_new+solution[n+step-1]);
    solution.push_back(U_new);
}

void Adams_Bashforth::set_beta(int p){
    if(p==1){
        beta[step-1]=1.0;
    }
    else if(p==2){
        beta[step-1]=1.5;
        beta[step-2]=-0.5;
    }
    else if(p==3){
        beta[step-1]=23.0/12.0;
        beta[step-2]=-16.0/12.0;
        beta[step-3]=5.0/12.0;
    }
    else if(p==4){
        beta[step-1]=55.0/24.0;
        beta[step-2]=-59.0/24.0;
        beta[step-3]=37.0/24.0;
        beta[step-4]=-9.0/24.0;
    }
    else{
        cerr<<"can't calculate!"<<endl;
    }
}

Adams_Moulton::Adams_Moulton(){
    method=Method::AMF;
    this->set_beta(1);
}

Adams_Moulton::Adams_Moulton(int p){
    method =Method::AMF;
    accuracy=p;
    this->set_beta(p);
}

void Adams_Moulton::set_beta(int p){
    if(p==2){
        beta[step]=0.5;
        beta[step-1]=0.5;
    }
    else if(p==3){
        beta[step]=5.0/12.0;
        beta[step-1]=8.0/12.0;
        beta[step-2]=-1.0/12.0;
    }
    else if(p==4){
        beta[step]=9.0/24.0;
        beta[step-1]=19.0/24.0;
        beta[step-2]=-5.0/24.0;
        beta[step-3]=1.0/24/0;
    }
    else if(p==5){
        beta[step]=251.0/720.0;
        beta[step-1]=646.0/720.0;
        beta[step-2]=-264.0/720.0;
        beta[step-3]=106.0/720.0;
        beta[step-4]=-19.0/720.0;
    }
}

void Adams_Moulton::OneStep(const Function &f, int n){
    Matrix U_new(1, d);
    Matrix U_old(1, d);
    U_old=solution[n+step-1];
    while (abs(U_new.l2_norm()-U_old.l2_norm())>1e-6){
        Matrix temp;
        temp=U_old;
        U_old=U_new;
        U_new=move(temp);
        U_new=move(U_new+f(U_new, (n+step)*k)*beta[n+step]);
        for(int j=0; j<step; ++j){
            U_new=move(U_new+f(solution[n+j], (n+j)*k)*(beta[j]*k+t_begin));
        }
    }
    solution.push_back(U_new);    
}

Backward_Differential::Backward_Differential(){
    method=Method::BDF;
    accuracy=1;
    this->set(1);
}

Backward_Differential::Backward_Differential(int p){
    method=Method::BDF;
    accuracy=p;
    this->set(p);
}

void Backward_Differential::set(int p){
    if(p==1){
        alpha[step-1]=-1.0;
        beta[step]=1.0;
    }
    else if(p==2){
        alpha[step-1]=-4.0/3.0;
        alpha[step-2]=1.0/3.0;
        beta[step]=2.0/3.0;
    }
    else if(p==3){
        alpha[step-1]=-18.0/11.0;
        alpha[step-2]=9.0/11.0;
        alpha[step-3]=-2.0/11.0;
        beta[step]=6.0/11.0;
    }
    else if(p==4){
        alpha[step-1]=-48.0/25.0;
        alpha[step-2]=36.0/25.0;
        alpha[step-3]=-16.0/25.0;
        alpha[step-4]=3.0/25.0;
        beta[step]=12.0/25.0;
    }
}

void Backward_Differential::OneStep(const Function &f, int n){
    Matrix U_new(1, d);
    Matrix U_old(1, d);
    U_old=solution[n+step-1];
    while (abs(U_new.l2_norm()-U_old.l2_norm())>1e-6){
        Matrix temp;
        temp=U_old;
        U_old=U_new;
        U_new=move(temp);
        U_new=move(U_new+f(U_new, (n+step)*k)*(k*beta[step]+t_begin));
        for(int j=0; j<step; ++j){
            U_new=move(U_new+solution[n+j]*(-alpha[j]));
        }
    }
    solution.push_back(U_new);    
}