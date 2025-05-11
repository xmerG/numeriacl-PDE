#include"Runge_Kutta.h"

void Implicit_Runge_Kutta::solve(const Function &f){
    for(int i=0; i<maxStep; ++i){
        this->OneStep(f, i);
    }
}

void Implicit_Runge_Kutta::OneStep(const Function &f, int n){
    Matrix y(s, d);
    for(int i=0; i<100; ++i){
        
    }
}


void classical_RK::OneStep(const Function &f, int n){
    double current_t=n*k+t_begin;
    Matrix y1=move(f(solution[n], current_t));
    current_t+=k/2.0;
    Matrix y2=move(f(solution[n]+y1*(k/2.0), current_t));
    Matrix y3=move(f(solution[n]+y2*(k/2.0), current_t));
    current_t+=k/2.0;
    Matrix y4=move(f(solution[n]+y3*k, current_t));
    Matrix U_new=move(solution[n]+(y1+y2*2+y3*2+y4)*(k/6.0));
    solution.push_back(U_new);
}


