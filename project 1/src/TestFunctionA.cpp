#include"TestFunctionA.hpp"

const double b=cos(1.0);
const double a=sin(1.0);

double  Laplacian::operator()(double x, double y) const{
    return (-1.0+sin(x)-pow(cos(x), 2))*exp(y+sin(x));
}
    






double DirichletF::operator()(double x, double y) const{
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
   
double NeumannF::operator()(double x, double y) const{
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
    else {
        cerr<<"not defined"<<endl;
        return -1;
    }
}

Mixed::Mixed(const vector<int> &_v):mixed(_v){}

double Mixed::operator()(double x, double y) const{
    if(x==0.0){
        if(mixed[1]==0){return exp(y);}
        else{return -exp(y);}
         
    }
    else if(x==1.0){
        if(mixed[2]==0){return exp(y+a);}
        else{return b*exp(y+a);}
    }
    else if(y==0.0){
        if(mixed[0]==0){return exp(sin(x));}
        else{return -exp(sin(x));}        
    }
    else if(y==1.0){
        if(mixed[3]==0){return exp(1.0+sin(x));}
        else{return exp(sin(x)+1);}
    }
    else {
        cerr<<"not defined"<<endl;
        return -1;
    }
}
    
double primitive::operator()(double x, double y) const{
    return exp(y+sin(x));
}

irreNeumann::irreNeumann(){}

irreNeumann::irreNeumann(Circle *_c):c(_c){}

double irreNeumann::operator()(double x, double y) const{
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

irreMixed::irreMixed(){}

irreMixed::irreMixed(Circle *_c, const vector<int> &_v):c(_c), v(_v){}

double irreMixed::operator()(double x, double y) const{
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
    else if(c->onCircle(x,y)){
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


    