#ifndef _MATRIX_
#define _MATRIX_
#include<vector>
#include<iostream>
using namespace std;

class IVP_solver;
class Function;

class Matrix{
private:
    int m=0; //rows
    int n=0; //columns
    vector<vector<double>> elements;
public:
    Matrix();
    Matrix(const int &, const int &);
    Matrix(const int &, const int &, const vector<vector<double>> &);
    Matrix(const vector<vector<double>> &);
    Matrix(Matrix &&other);
    Matrix operator=(Matrix &&other);
    double operator()(const int &i, const int&j);
    Matrix operator*(const Matrix &other) const;
    Matrix operator+(const Matrix &other) const;
    Matrix operator-(const Matrix &other) const;
    void add_elements(const Matrix &other);
    void set_elements(const Matrix &,const int &);
    void print() const;
    friend IVP_solver;
    friend Function;
};

#endif