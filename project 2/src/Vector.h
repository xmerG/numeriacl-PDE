#ifndef _VECTOR_
#define _VECTOR_
#include<vector>
#include<iostream>
using namespace std;

class Sparse_Matrix;

class Vector{
private:
    int n=0;
    vector<double> elements;
public:
    Vector();
    Vector(const int &n);
    Vector(const int &n, const vector<double> &_e);
    Vector operator+(const Vector &v);
    Vector operator-(const Vector &v);
    Vector operator*(const double &a);
    friend class Sparse_Matrix;
};

#endif
