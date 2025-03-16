#ifndef TESTFUNCTIONA
#define TESTFUNCTIONA
#include"Circle.hpp"
#include<vector>
#include"Function.hpp"
#include"BoundaryCondition.hpp"
#include"Domain.hpp"
using namespace std;

/*
class boundaryFunction:public Function{
private:
    Circle *c;
    BoundaryCondition bc;
    Domain d;
    vector<int> v;
public:
    boundaryFunction(BoundaryCondition _bc, Domain _d, Circle *_c=nullptr,const vector<int> &_v=vector<int>{});
    double operator()(double x, double y) const;
};*/



class Laplacian:public Function{
public:
    double operator()(double x, double y) const;
};

class DirichletF:public Function{
public:
    double operator()(double x, double y) const;
};

class NeumannF:public Function{
public:
    double operator()(double x, double y) const;
};

class Mixed:public Function{
private:
    vector<int> mixed;
public:
    Mixed(const vector<int> &_v);
    double operator()(double x, double y) const;
};


class primitive:public Function{
public:
    double operator()(double x, double y) const;
};

class irreMixed:public Function{
private:
    Circle *c=nullptr;
    vector<int> v;
public:
    irreMixed();
    irreMixed(Circle *_c,const vector<int> &_v);
double operator()(double x, double y) const;
};

class irreNeumann:public Function{
private: 
    Circle *c=nullptr;
public:
    irreNeumann();
    irreNeumann(Circle *_c);
    double operator()(double x, double y) const;
};

#endif