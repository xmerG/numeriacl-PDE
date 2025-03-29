#include"Vector.h"

Vector::Vector(){}

Vector::Vector(const int &_n):n(_n){
    elements.resize(n, 0.0);
    m=sqrt(n);
}

Vector::Vector(const vector<double> &e):elements(e){n=e.size(); m=sqrt(n);}

Vector::Vector(const int &_n, const vector<double> & _e):n(_n), elements(_e){m=sqrt(n);}

Vector::Vector(Vector &&other) noexcept: n(other.n), m(other.m),elements(std::move(other.elements)) {
    other.n = 0;
    other.m=0;
}

Vector& Vector::operator=(Vector&& other) noexcept {
    if (this != &other) {
        n = other.n;
        m=other.m;
        elements = std::move(other.elements);
        other.n = 0;
        other.m=0;
    }
    return *this;
}

void Vector::set_Value(const int &i, const double &value){
    if(i>=0 && i<n){
        this->elements[i]=value;
    }
}

Vector Vector::operator+(const Vector &v) const{
    if(n!=v.n){
        cerr<<"invalid!"<<endl;
        return Vector();
    }
    else{
        vector<double> new_elements(n, 0.0);
        for(int i=0; i<n; ++i){
            new_elements[i]=elements[i]+v.elements[i];
        }
        return Vector(n, new_elements);
    }
}

Vector Vector::operator-(const Vector &v) const{
    if(n!=v.n){
        cerr<<"invalid!"<<endl;
        return Vector();
    }
    else{
        vector<double> new_elements(n, 0.0);
        for(int i=0; i<n; ++i){
            new_elements[i]=elements[i]-v.elements[i];
        }
        return Vector(n, new_elements);
    }
}

Vector Vector::operator*(const double &a) const{
    vector<double> new_elements(n,0.0);
    for(int i=0; i<n; ++i){
        new_elements[i]=a*elements[i];
    }
    return Vector(n, new_elements);
}


double Vector::operator()(const int &i) const{
    if(i>=0 && i<n){
        return this->elements[i];
    }
    else{
        return 0.0;
    }
}

int Vector::getdim() const{return this->n;}

void Vector::print() const{
    for(int i=0; i<n; ++i){
        cout<<elements[i]<<" ";
    }
    cout<<endl;
}

double Vector::operator()(const int &i, const int &j) const{
    if(i>=0 && i<m && j>=0 && j<m){
        return elements[i+j*m];
    }
    else{return 0.0;}
}

void Vector::set_Value(const int &i, const int &j, const double &value){
    if(i>=0 && i<m && j>=0 && j<m){
        elements[i+j*m]=value;
    }
}