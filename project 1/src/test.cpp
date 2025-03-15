#include<iostream>
#include"EquationSolver.hpp"
#include<cmath>
#include"Circle.hpp"
using namespace std;

const double b=cos(1.0);
const double a=sin(1.0);
class Laplacian:public Function{
public:
    double operator()(double x, double y) const{
        return (-1.0+sin(x)-pow(cos(x), 2))*exp(y+sin(x));
    }
};

class DirichletF:public Function{
public:
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
public:
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

class Mixed:public Function{
public:
    double operator()(double x, double y) const{
        if(x==0.0){
            return exp(y); 
        }
        else if(x==1.0){
            return exp(y+a);
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
public:
    double operator()(double x, double y) const{
        return exp(y+sin(x));
    }
};

class irreNeumann:public Function{
private: 
    Circle *c=nullptr;
public:
    irreNeumann(){}
    irreNeumann(Circle *_c):c(_c){}
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
        else if(c->onCircle(x,y)){
            double r=c->get_radius();
            double x0=c->getX();
            double y0=c->getY();
            return ((x0-x)*cos(x)+(y0-y))*exp(y+sin(x))/r;
        }
        else {
            cerr<<"not defined"<<endl;
            return -1;
        }
    }
};

int main(){
    Laplacian f1;
    DirichletF g;
    primitive f0;

    /*EquationSolver<Domain::regular, BoundaryCondition::Dirichlet> solver1(32, f1);
    solver1.solveEquation(g);
    solver1.norm_error(f0);
    solver1.print("test.json", f0);*/

    /*vector<int> neumann=vector<int>{1,1,1,1};
    NeumannF g1;
    EquationSolver<Domain::regular, BoundaryCondition::Neumann> solver2(4,f1);
    solver2.solveEquation(g1,f0(1.0/4.0, 1.0/4.0),neumann);
    solver2.norm_error(f0);
    solver2.print("test.json", f0);*/

    /*vector<int> mixed=vector<int>{1,0,0,1};
    Mixed h1;
    EquationSolver<Domain::regular, BoundaryCondition::Mixed> solver5(32, f1);
    solver5.solveEquation(h1, 0.0,mixed);
    solver5.norm_error(f0);
    solver5.print("test.json", f0);*/



    vector<double> D{f0(0.0,0.0), f0(1.0,0.0), f0(0.0,1.0), f0(1.0,1.0)};
    Circle *c = new Circle(0.5, 0.5, 0.2);
    /*EquationSolver<Domain::irregular,BoundaryCondition::Dirichlet> solver3(8,f1,c);
    solver3.solveEquation(f0);
    solver3.print("test.json", f0);
    delete c;*/

    irreNeumann g1(c);
    EquationSolver<Domain::irregular, BoundaryCondition::Neumann> solver4(64,f1,c);
    solver4.solveEquation(g1, f0(1.0/64.0,1.0/64.0));
    solver4.print("test.json", f0);
    solver4.norm_error(f0);
    return 0;
}