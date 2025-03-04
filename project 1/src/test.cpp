#include<iostream>
#include"EquationSolver.hpp"
#include<cmath>
using namespace std;

const double b=cos(1.0);
const double a=sin(1.0);
class Laplacian:public Function{
    double operator()(double x, double y) const{
        return (-1.0+sin(x)-pow(cos(x), 2))*exp(y+sin(x));
    }
};

class DirichletF:public Function{
    double operator()(double x, double y) const{
        if(x==0.0){
            return exp(y); 
        }
        else if(x==1.0){
            return exp(y+a);
        }
        else if(y==0.0){
            return exp(sin(x));
        }
        else if(y==1.0){
            return exp(1.0+sin(x));
        }
        else{
            cerr<<"not defined"<<endl;
            return -1;
        }
    }
};

class NeumannF:public Function{
    double operator()(double x, double y) const{
        if(x==0.0){
            return -exp(y);
        }
        else if(x==1.0){
            return b*exp(y+a);
        }
        else if(y==0.0){
            return -exp(sin(x));
        }
        else if(y==1.0){
            return exp(sin(x)+1);
        }
        else {
            cerr<<"not defined"<<endl;
            return -1;
        }
    }
};

class primitive:public Function{
    double operator()(double x, double y) const{
        return exp(y+sin(x));
    }
};

int main(){
    Laplacian f1;
    //DirichletF g;
    primitive f0;
    //EquationSolver<Domain::regular, BoundaryCondition::Dirichlet> solver1(4, f1);
    //solver1.solveEquation(g);
    //solver1.norm_error(f0);
    //solver1.print("test.json", f0);

    NeumannF g1;
    EquationSolver<Domain::regular, BoundaryCondition::Neumann> solver2(4,f1);
    solver2.solveEquation(g1);
    solver2.print("test.json", f0);
    return 0;
}