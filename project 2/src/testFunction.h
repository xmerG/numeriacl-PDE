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

class Neumann1:public Function{
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

class Mixed1:public Function{
private:
    vector<int> mixed;
public:
    Mixed1(){}
    Mixed1(const vector<int> &v):mixed(v){}
    double operator()(const double &x, const double &y) const {return 0.0;}
    //--------------------------------------------------------
    double operator()(const double &x)const{
        if(x==0){
            if(mixed[0]==0){
                return 1.0;
            }
            else{
                return -2.0;
            }
        }

        else if(x==1.0){
            if(mixed[1]==0){
                return exp(1.0+a);
            }
            else{
                return (1.0+b)*exp(1.0+a);
            }
        }
        else{
            cerr<<"not defined"<<endl;
            return 0.0;
        }
    }
};

class Laplacian1:public Function{
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
        //直接给出角上的函数值
        if(x==0.0 && y==0.0){
            //return 1.0;
            return -2.0;
        }
        else if(x==0.0 && y==1.0){
            //return exp(1.0);
            return 0.0;
        }
        else if(x==1.0 && y==0.0){
            //return exp(a);
            return (b-1.0)*exp(a);
        }
        else if(x==1.0 && y==1.0){
            //return exp(a+1.0);
            return (b+1.0)*exp(a+1.0);
        }


        else{
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

    }
    double operator()(const double &x)const{return 0.0;}
};

class Mixed2:public Function{
private:
    vector<int> mixed=vector<int>{0,0,0,0};
public:
    Mixed2(){}
    Mixed2(const vector<int> &_v):mixed(_v){}
    double operator()(const double &x)const{return 0.0;}
    double operator()(const double &x, const double &y) const {
        //直接给出角上的函数值
        if(x==0.0 && y==0.0){
            return 1.0;
        }
        else if(x==0.0 && y==1.0){
            return d;
        }
        else if(x==1.0 && y==0.0){
            return exp(a);
        }
        else if(x==1.0 && y==1.0){
            return exp(a+1.0);
        }


        else{
            if(x==0.0){
                if(mixed[1]==1){
                    return -exp(y);
                }
                else{
                    return exp(y);
                }
            }
            else if(x==1.0){
                if(mixed[2]==1){
                    return b*exp(y+a);
                }
                else{
                    return exp(y+a);
                }
            }
            else if(y==0.0){
                if(mixed[0]==1){
                    return -exp(sin(x));
                }
                else{
                    return exp(sin(x));
                }
            }
            else if(y==1.0){
                return exp(sin(x)+1.0);
            }
            else{
                cerr<<"not defined on interior"<<endl;
                return 0.0;
            }
        }

    }
};

#endif