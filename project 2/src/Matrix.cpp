#include"Matrix.hpp"
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

//void Sparse_Matrix::transform(){
    //for(map<label, double>::iterator i=A.begin(); i!=A.end(); ++i){
        
    //}
//}

Sparse_Matrix Sparse_Matrix::operator+(const Sparse_Matrix &B){
    if(n!=B.n){
        cerr<<"invalid since the matrix dimension must agree for addision"<<endl;
        return;
    } 
    else{
        map<label, double> C;
        map<label, double> A1=B.A;
        for(map<label, double>::iterator i=A1.begin(); i!=A1.end(); ++i){
            if(A.count(i->first)){
                C[i->first]=i->second+A[i->first];
            }
            else{
                C[i->first]=i->second;
            }
        }
        return Sparse_Matrix(this->n, C);
    }
}

Sparse_Matrix Sparse_Matrix::operator*(const double &a){
    map<label, double> elements;
    for(map<label,double>::iterator i=A.begin(); i!=A.end(); ++i){
        elements[i->first]=i->second*a;
    }
    return Sparse_Matrix(this->n, elements);
}

Vector Sparse_Matrix::operator*(const Vector &v){
    vector<double> elements;
    

}