#ifndef _MATRIX_
#define _MATRIX_
#include<vector>
#include<iostream>
#include<cmath>
using namespace std;

class IVP_solver;
class Function;

class Matrix{
protected:
    int m=0; //rows
    int n=0; //columns
    vector<vector<double>> elements;
public:
    Matrix();
    Matrix(const int &, const int &);
    Matrix(const int &, const int &, const vector<vector<double>> &);
    Matrix(const vector<vector<double>> &);
    Matrix(Matrix &&other);
    Matrix operator=(const Matrix& other);
    Matrix operator=(Matrix &&other);
    int get_row() const;
    int get_col() const;
    double operator()(const int &i, const int&j);
    Matrix operator()(const int &, const int &, const int &, const int &);
    Matrix operator*(const double & );
    Matrix operator*(const Matrix &other) const;
    Matrix operator+(const Matrix &other) const;
    Matrix operator-(const Matrix &other) const;
    void add_row_elements(const Matrix &other);
    void add_col_elements(const Matrix &other);
    void set_elements(const int &, const int &, const double &);
    void set_elements(const Matrix &, const int &, const int &, const int &, const int&);
    vector<Matrix> LU() const;
    Matrix operator/(Matrix &);
    double l2_norm() const;
    void print() const;
    friend IVP_solver;
    friend Function;
};

#endif