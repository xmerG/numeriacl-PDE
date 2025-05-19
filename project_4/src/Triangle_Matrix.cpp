#include"Triangle_Matrix.h"
#include<cmath>

Triangle_Matrix::Triangle_Matrix(const int &_n):Sparse_Matrix(n){}

Triangle_Matrix::Triangle_Matrix(Triangle_Matrix &&other){
    n=other.n;
    elements=move(other.elements);
}

Triangle_Matrix Triangle_Matrix::Cholesky() {
    Triangle_Matrix L(n);
    L.setValues(0, 0, sqrt(this->operator()(0,0)));
    L.setValues(1, 0, elements[1][0]/L.elements[0][0]);
    for(int i=1; i<this->n; ++i){
        double value=sqrt(elements[i][i]-pow(L.elements[i][i-1], 2));
        L.setValues(i, i, value);
        value=elements[i+1][i]/L.elements[i][i];
        L.setValues(i+1, i, value);
    }
    return move(L);
}

Triangle_Matrix& Triangle_Matrix::operator=(Sparse_Matrix &&other){
    Sparse_Matrix::operator=(other);
    return *this;
}

Vector Triangle_Matrix::operator/(const Vector &b){
    Triangle_Matrix L=this->Cholesky();
    Vector sol(n);
    vector<double> b_elements=b.get_elements();
    sol.set_value(n-1, b_elements[n-1]/L.elements[n-1][n-1]);
    for(int i=n-2; i>=0; ++i){
        double value=(b_elements[i]-L.elements[i+1][i]*sol(i+1))/L.elements[i][i];
        sol.set_value(i, value);
    }
    sol.set_value(0, L.elements[0][0]);
    for(int i=1; i<n; ++i){
        double value=(sol(i)-sol(i+1)*L.elements[i][i-1])/L.elements[i][i];
        sol.set_value(i, value);
    }
    return sol;
}

Triangle_Matrix& Triangle_Matrix::operator=(Triangle_Matrix &&other){
    if (this != &other) {
        n = other.n;
        elements = move(other.elements);
        other.n = 0;
    }
    return *this;
}
