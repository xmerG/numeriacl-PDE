#ifndef _TESTFUNCTION_
#define _TESTFUNCTION_
#include"Function.hpp"
#include<cmath>

const double b=cos(1.0);
const double a=sin(1.0);
const double d=exp(1.0);

class F1:public Function{
public:
    double operator() (const double &x, const double &y) const{return 0.0;}
    double operator()(const double &x) const{
        return exp(x+sin(x));
    }
};

class Neumann:public Function{
public:
    double operator()(const double &x, const double &y) const {return 0.0;}
    double operator()(const double &x)const{
        if(x==0){
            return -2.0;
        }
        else if(x==1.0){
            return (1.0+cos(1.0))*exp(1.0+sin(1.0));
        }
        else{return -1.0;}
    }
};

class Laplacian:public Function{
public:
    double operator()(const double &x, const double &y) const {return 0.0;}
    double operator()(const double &x)const{
        return exp(x+sin(x))*(sin(x)-(1+cos(x))*(1+cos(x)));
    }
};


class F2:public Function{
public:
    double operator() (const double &x, const double &y) const{
        return exp(y+sin(x));
    }
    double operator()(const double &x) const{return 0.0;}
};

class Laplacian2:public Function{
public:
    double operator()(const double &x, const double &y) const {
        return (-1.0+sin(x)-pow(cos(x), 2))*exp(y+sin(x));
    }
    double operator()(const double &x)const{return 0.0;}
};

class Neumann2:public Function{
public:
    double operator()(const double &x, const double &y) const {
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
            return exp(sin(x)+1.0);
        }
        else{
            cerr<<"not defined on interior"<<endl;
            return 0.0;
        }
    }
    double operator()(const double &x)const{return 0.0;}
};

#endif