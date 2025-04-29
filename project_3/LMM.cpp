#include"LMM.h"

void LMM::solve(const Function &f){
    this->setInitials(f);
    for(int i=0; i<maxStep-step; ++i){
        this->OneStep(f, i);
    }
}

void LMM::solve(const Function &f, const double &_initial, const double &t1, const double &t2, int N){
    IVP_solver::setValues(t1, t2, N);
    LMM::solve(f);
}

void LMM::setInitials(const Function &f){
    
}