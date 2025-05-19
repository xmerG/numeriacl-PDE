#include"Vector.h"

Vector::Vector(const int &_n):n(_n){
    elements.resize(n, 0.0);
}

Vector::Vector(const vector<double> &_elements):elements(_elements){
    n=elements.size();
}

double Vector::operator()(const int &index){
    if(index>=0 && index<n){
        return elements[index];
    }
    else{
        cerr<<"index out of range!"<<endl;
        return 0.0;
    }
}

Vector& Vector::operator=(const Vector &other) {
    if (this != &other) {  
        n = other.n;
        elements = other.elements;  
    }
    return *this;  
}

Vector& Vector::operator=(Vector &&other) noexcept {
    if (this != &other) {  
        n = std::move(other.n);
        elements = std::move(other.elements);  
        other.n = 0;
    }
    return *this;
}

void Vector::set_value(const int &i, const double &value){
    elements[i]=value;
}

Vector::Vector(const Vector &other):n(other.n), elements(other.elements){}

Vector::Vector(Vector &&other):n(other.n), elements(move(other.elements)){}

Vector Vector::operator+(const Vector &other){
    Vector v(this->n);
    if(this->n == other.n){
        for(int i=0; i<n; ++i){
            v.elements[i]=this->elements[i]+other.elements[i];
        }
    }
    else{
        cerr<<"dimensions don't match!"<<endl;
    }
    return move(v);
}

Vector Vector::operator-(const Vector &other){
    Vector v(this->n);
    if(this->n == other.n){
        for(int i=0; i<n; ++i){
            v.elements[i]=this->elements[i]-other.elements[i];
        }
    }
    else{
        cerr<<"dimensions don't match!"<<endl;
    }
    return move(v);    
}

Vector Vector::operator*(const double &a){
    Vector v(n);
    for(int i=0; i<n; ++i){
        v.elements[i]=this->elements[i]*a;
    }
    return move(v);
}

vector<double> Vector::get_elements() const{
    return this->elements;
}