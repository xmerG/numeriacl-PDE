#ifndef _MATRIX_
#define _MATRIX_
#include<utility>
#include<map>
#include"Vector.h"
using namespace std;
using label=pair<int, int>;

//only valid for square matrix
class Sparse_Matrix{
private: 
    map<label, double> A;
    int n;
    
public:
    Sparse_Matrix();
    Sparse_Matrix(const int &_n);
    Sparse_Matrix(const int &_n, const map<label, double> &A);
    void setValues(const int &i, const int &j, const double &value); //set A(i,j)
    double operator()(const int &i, const int &j);
    void transform();
    Sparse_Matrix operator+(const Sparse_Matrix &B);
    Sparse_Matrix operator*(const Sparse_Matrix &B);
    Sparse_Matrix operator*(const double &a);
    Vector operator*(const Vector &v);

};

#endif