#include"Sparse_Matrix.h"
Sparse_Matrix::Sparse_Matrix():n(0){}

Sparse_Matrix::Sparse_Matrix(const int &_n):n(_n){}

Sparse_Matrix::Sparse_Matrix(const int &_n, const map<label, double> &_A):n(_n), A(_A){}

void Sparse_Matrix::setValues(const int &i, const int &j, const double &value){
    label l(i,j);
    if(i>=0 && i<n && j>=0 && j<n){
        A[l]=value;
    }
    else{
        cerr<<"invalid label!"<<endl;
    }
}

double Sparse_Matrix::operator()(const int &i, const int &j){
    if(i>=0 && i<n && j>=0 && j<n){
        label l(i,j);
        return A[l];
    }
}

Sparse_Matrix Sparse_Matrix::operator+(const Sparse_Matrix &B){
    if(n!=B.n){
        cerr<<"invalid since the matrix dimension must agree for addision"<<endl;
    } 
    //-------------------不能直接修改A-------------------------
    else{
        Sparse_Matrix C(n);
        map<label, double> A1=B.A;
        for(map<label, double>::iterator i=A1.begin(); i!=A1.end(); ++i){
            if(A.count(i->first)){
                A[i->first]+=i->second;
            }
            else{
                A[i->first]=i->second;
            }
        }
    }
}
