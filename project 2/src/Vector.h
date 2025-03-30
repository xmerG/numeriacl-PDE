#ifndef _VECTOR_
#define _VECTOR_
#include<vector>
#include<iostream>
#include<cmath>
#include <algorithm>
#include<fstream>
#include<string>
#include <nlohmann/json.hpp>
using namespace std;

class Sparse_Matrix;

class Vector{
private:
    int n=0;
    int m=0;
    vector<double> elements;
public:
    Vector();
    Vector(const int &n);
    Vector(const vector<double> &e);
    Vector(const int &n, const vector<double> &_e);
    Vector(Vector&& other) noexcept;
    Vector(const Vector&) = delete;           // 禁用拷贝
    Vector& operator=(const Vector&) = delete;
    void set_Value(const int &i, const double &value);
    void set_Value(const int &i, const int &j, const double &value);
    Vector& operator=(Vector&& other) noexcept;
    double operator()(const int &i) const;
    Vector operator+(const Vector &v) const;
    Vector operator-(const Vector &v) const;
    Vector operator*(const double &a) const;
    int getdim() const;
    void print() const;
    double operator()(const int &i, const int &j) const;
    void go_zero();
    friend class Sparse_Matrix;
    vector<double> getelements() const;
    void print_to_file(const string &filename);
};

#endif
