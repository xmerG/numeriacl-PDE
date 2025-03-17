#include"TestFunction.hpp"
using namespace std;
const double b=cos(1.0);
const double a=sin(1.0);
const double d=exp(1.0);

//TEST A
double LaplacianA::operator()(double x,double y) const{
    return (-1.0+sin(x)-pow(cos(x), 2))*exp(y+sin(x));
}


double PrimitiveA::operator()(double x, double y) const{
    return exp(y+sin(x));
}

NeumannA::NeumannA(){}
NeumannA::NeumannA(Circle *_c):c(_c){}
double NeumannA::operator()(double x, double y) const{
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
    else if(c!=nullptr && c->onCircle(x,y)){
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

MixedA::MixedA(){}
MixedA::MixedA(Circle *_c, const vector<int> &_v):c(_c), v(_v){}
double MixedA::operator()(double x, double y) const{
    if(x==0.0){
        if(v[1]==0){return exp(y);}
        else{return -exp(y);}
         
    }
    else if(x==1.0){
        if(v[2]==0){return exp(y+a);}
        else{return b*exp(y+a);}
    }
    else if(y==0.0){
        if(v[0]==0){return exp(sin(x));}
        else{return -exp(sin(x));}        
    }
    else if(y==1.0){
        if(v[3]==0){return exp(1.0+sin(x));}
        else{return exp(sin(x)+1);}
    }
    else if(c!=nullptr && c->onCircle(x,y)){
        if(v[4]==1){
            double r=c->get_radius();
            double x0=c->getX();
            double y0=c->getY();
            return ((x0-x)*cos(x)+(y0-y))*exp(y+sin(x))/r;
        }
        else{
            return exp(y+sin(x));
        }
    }
    else {
        cerr<<"not defined"<<endl;
        return -1;
    }
}

//TEST B
double LaplacianB::operator()(double x, double y) const{
    return 5.0/pow(x+2*y+0.5, 2);
}

double PrimitiveB::operator()(double x, double y) const{
    return log(0.5+x+2*y);
}

NeumannB::NeumannB(){}
NeumannB::NeumannB(Circle *_c):c(_c){}

double NeumannB::operator()(double x, double y) const
{
    if(x==0.0){
        return -1.0/(2*y+0.5);
    }
    else if(x==1.0){
        return 1.0/(2*y+1.5);
    }
    else if(y==0.0){
        return -2.0/(0.5+x);
    }
    else if(y==1.0){
        return 2.0/(2.5+x);
    }
    else if(c!=nullptr && c->onCircle(x,y)){
        double r=c->get_radius();
        double x0=c->getX();
        double y0=c->getY();
        return ((x0-x)+2*(y0-y))/(r*(x+2*y+0.5));
    }
    else {
        cerr<<"not defined"<<endl;
        return -1;
    }
}

MixedB::MixedB(){}
MixedB::MixedB(Circle *_c, const vector<int> &_v):c(_c), v(_v){}
double MixedB::operator()(double x, double y)const{
    if(x==0.0){
        if(v[1]==0){return log(0.5+2*y);}
        else{return -1.0/(2*y+0.5);}
         
    }
    else if(x==1.0){
        if(v[2]==0){return log(1.5+2*y);}
        else{return 1.0/(2*y+1.5);}
    }
    else if(y==0.0){
        if(v[0]==0){return log(0.5+x);}
        else{return -2.0/(0.5+x);}        
    }
    else if(y==1.0){
        if(v[3]==0){return log(2.5+x);}
        else{return 2.0/(2.5+x);}
    }
    else if(c!=nullptr && c->onCircle(x, y)){
        if(v[4]==0){return log(0.5+x+2*y);}
        else{
            double r=c->get_radius();
            double x0=c->getX();
            double y0=c->getY();
            return ((x0-x)+2*(y0-y))/(r*(x+2*y+0.5));
        }
    }
    else {
        cerr<<"not defined"<<endl;
        return -1;
    }
}


//TEST C
double LaplacianC::operator()(double x, double y) const{
    double temp=exp(x)+y+1;
    return sin(temp)*(exp(2*x)+1)-cos(temp)*exp(x);
}

double PrimitiveC::operator()(double x, double y) const{
    return sin(exp(x)+y+1);
}

NeumannC::NeumannC(){}
NeumannC::NeumannC(Circle *_c):c(_c){}
double NeumannC::operator()(double x, double y) const{
    if(x==0.0){
        return -cos(y+2.0);
    }
    else if(x==1.0){
        return cos(d+1.0+y)*d;
    }
    else if(y==0.0){
        return -cos(exp(x)+1.0);
    }
    else if(y==1.0){
        return cos(exp(x)+2.0);
    }
    else if(c!=nullptr && c->onCircle(x,y)){
        double r=c->get_radius();
        double x0=c->getX();
        double y0=c->getY();
        return ((x0-x)*exp(x)+(y0-y))*cos(exp(x)+y+1)/r;
    }
    else {
        cerr<<"not defined"<<endl;
        return -1;
    }
}

MixedC::MixedC(){}
MixedC::MixedC(Circle *_c, const vector<int> &_v):c(_c), v(_v){}
double MixedC::operator()(double x, double y)const{
    if(x==0.0){
        if(v[1]==0){return sin(y+2.0);}
        else{return -cos(y+2.0);}
         
    }
    else if(x==1.0){
        if(v[2]==0){return sin(d+y+1.0);}
        else{return cos(d+1+y)*d;}
    }
    else if(y==0.0){
        if(v[0]==0){return sin(exp(x)+1.0);}
        else{return -cos(exp(x)+1);}        
    }
    else if(y==1.0){
        if(v[3]==0){return sin(exp(x)+2.0);}
        else{return cos(exp(x)+2.0);}
    }
    else if(c!=nullptr && c->onCircle(x, y)){
        if(v[4]==0){return sin(exp(x)+y+1);}
        else{
            double r=c->get_radius();
            double x0=c->getX();
            double y0=c->getY();
            return ((x0-x)*exp(x)+(y0-y))*cos(exp(x)+y+1)/r;
        }
    }
    else {
        cerr<<"not defined"<<endl;
        return -1;
    }
}


