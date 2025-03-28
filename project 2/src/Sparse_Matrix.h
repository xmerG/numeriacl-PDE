#ifndef _SPARSE_MATRIX_
#define _SPARSE_MATRIX_
#include<utility>
#include<map>
#include"Vector.h"
using namespace std;
using label=map<int, double>; //the first int denote jth column while the second int denote the element

//only valid for square matrix
class Sparse_Matrix{
private: 
    vector<label> elements;
    int n;
public:
    Sparse_Matrix();
    Sparse_Matrix(const int &_n);
    Sparse_Matrix(const int &_n, const vector<label> &e);
    Sparse_Matrix(Sparse_Matrix&& other) noexcept;
    Sparse_Matrix& operator=(Sparse_Matrix&& other) noexcept;
    void setValues(const int &i, const int &j, const double &value); //set A(i,j)
    double operator()(const int &i, const int &j) const;
    void transform();
    Sparse_Matrix operator+(const Sparse_Matrix &B) const;
    Sparse_Matrix operator*(const Sparse_Matrix &B) const;
    Sparse_Matrix operator*(const double &a) const;
    Vector operator*(const Vector &v) const;

};

#endif