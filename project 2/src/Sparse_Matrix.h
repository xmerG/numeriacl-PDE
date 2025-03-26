#ifndef _SPARSE_MATRIX_
#define _SPARSE_MATRIX_
#include<utility>
#include<iostream>
#include<map>
using namespace std;
using label=pair<int, int>;

class Sparse_Matrix{
private: 
    map<label, double> A;
    int n;
public:
    Sparse_Matrix();
    Sparse_Matrix(const int &_n);
    void setValues(const int &i, const int &j, const double &value); //set A(i,j)
    double operator()(const int &i, const int &j);
    Sparse_Matrix operator+()

};

#endif