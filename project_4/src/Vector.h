#ifndef _VECTOR_
#define _VECTOR_
#include<vector>
#include<iostream>
using namespace std;

class Vector{
private:
    vector<double> elements;
    int n;
public:
    Vector(const int &_n);
    Vector(const vector<double> &_elements);
    Vector(const Vector &other);
    Vector(Vector &&other);
    double operator()(const int &index);
    Vector& operator=(const Vector &other);
    Vector& operator=(Vector &&other);
    Vector operator*(const double & a);
    Vector operator-(const Vector &other);
    Vector operator+(const Vector &other);
    void set_value(const int &index, const double &value);
    vector<double> get_elements() const;
};

#endif