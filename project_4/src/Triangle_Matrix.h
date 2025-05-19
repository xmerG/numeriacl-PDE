#ifndef _TRIANGLE_MATRIX_
#define _TRIANGLE_MATRIX_
#include"Sparse_Matrix.h"
using namespace std;

class Triangle_Matrix:public Sparse_Matrix{
public:
    using Sparse_Matrix::operator=;
    Triangle_Matrix(Triangle_Matrix &&other);
    Triangle_Matrix(const int &n);
    Triangle_Matrix& operator=(Sparse_Matrix&& other);
    Triangle_Matrix& operator=(Triangle_Matrix && other);
    Triangle_Matrix Cholesky();
    Vector operator/(const Vector &v);
};

#endif