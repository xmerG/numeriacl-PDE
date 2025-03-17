#ifndef TESTFUNCTION
#define TESTFUNCTION
#include"Function.hpp"
#include<cmath>
#include<vector>
#include"Circle.hpp"
#include<iostream>


class LaplacianA : public Function {
public:
    double operator()(double x, double y) const;
};

class PrimitiveA:public Function{
public:
    double operator()(double x, double y) const;
};


class NeumannA:public Function{
private: 
    Circle *c=nullptr;
public:
    NeumannA();
    NeumannA(Circle *_c);
    double operator()(double x, double y) const;
};




class MixedA:public Function{
private:
    Circle *c=nullptr;
    vector<int> v;
public:
    MixedA();
    MixedA(Circle *_c,const vector<int> &_v);
    double operator()(double x, double y) const;
};




class LaplacianB : public Function {
public:
    double operator()(double x, double y) const;
};



class PrimitiveB:public Function{
public:
    double operator()(double x, double y) const;
};



class NeumannB:public Function{
private: 
    Circle *c=nullptr;
public:
    NeumannB();
    NeumannB(Circle *_c);
    double operator()(double x, double y) const;
};



class MixedB:public Function{
private:
    Circle *c=nullptr;
    vector<int> v;
public:
    MixedB();
    MixedB(Circle *_c,const vector<int> &_v);
    double operator()(double x, double y) const;
};


class LaplacianC : public Function {
public:
    double operator()(double x, double y) const;
};
    
class PrimitiveC:public Function{
public:
    double operator()(double x, double y) const;
};
    
    
class NeumannC:public Function{
private: 
    Circle *c=nullptr;
public:
    NeumannC();
    NeumannC(Circle *_c);
    double operator()(double x, double y) const;
};
    
    
    
    
class MixedC:public Function{
private:
    Circle *c=nullptr;
    vector<int> v;
public:
    MixedC();
    MixedC(Circle *_c,const vector<int> &_v);
    double operator()(double x, double y) const;
};

#endif