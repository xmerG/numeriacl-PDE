#include"Function.hpp"
#include<cmath>
class F1:public Function{
public:
    double operator() (const double &x, const double &y) const{return 0.0;}
    double operator()(const double &x) const{
        return exp(x+sin(x));
    }
};

class Laplacian:public Function{
public:
    double operator()(const double &x, const double &y) const {return 0.0;}
    double operator()(const double &x)const{
        return -(pow(1+cos(x), 2)-sin(x))*exp(x+sin(x));
    }
};