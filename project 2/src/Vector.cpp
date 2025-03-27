#include"Vector.h"

Vector::Vector(){}

Vector::Vector(const int &_n):n(_n){}

Vector::Vector(const int &_n, const vector<double> & _e):n(_n), elements(_e){}

Vector Vector::operator+(const Vector &v){
    if(n!=v.n){
        cerr<<"invalid!"<<endl;
    }
    else{
        vector<double> new_elements(n, 0.0);
        for(int i=0; i<n; ++i){
            new_elements[i]=elements[i]+v.elements[i];
        }
        return Vector(n, new_elements);
    }
}

Vector Vector::operator-(const Vector &v){
    if(n!=v.n){
        cerr<<"invalid!"<<endl;
    }
    else{
        vector<double> new_elements(n, 0.0);
        for(int i=0; i<n; ++i){
            new_elements[i]=elements[i]-v.elements[i];
        }
        return Vector(n, new_elements);
    }
}

Vector Vector::operator*(const double &a){
    vector<double> new_elements(n,0.0);
    for(int i=0; i<n; ++i){
        new_elements[i]=a*elements[i];
    }
    return Vector(n, new_elements);
}
