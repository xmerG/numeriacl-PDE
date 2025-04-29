#include"Matrix.h"

Matrix::Matrix(){}

Matrix::Matrix(const int &_m, const int &_n):m(_m), n(_n){
    elements.resize(m, vector<double>(n, 0.0));
}

Matrix::Matrix(const int &_m, const int &_n, const vector<vector<double>> &_elements):m(_m),
                n(_n), elements(_elements){}

Matrix::Matrix(const vector<vector<double>> &_elements):elements(_elements){
    m=elements.size();
    if(m!=0){
        n=elements[0].size();
    }
}

Matrix::Matrix(Matrix &&other):m(other.m), n(other.n),elements(move(other.elements)){}

Matrix Matrix::operator=(Matrix &&other){
    this->m=other.m;
    this->n=other.n;
    this->elements.clear();
    this->elements=move(other.elements);
    other.m=0;
    other.n=0;
}

double Matrix::operator()(const int &i, const int &j){
    if(0<=i<m && 0<=j<n){
        return elements[i][j];
    }
    else{
        return 0.0;
        cerr<<"index out of range!"<<endl;
    }
}

Matrix Matrix::operator*(const Matrix &other) const{
    Matrix newmatrix(this->m, other.n);
    if(this->n == other.m){
        for(int i=0; i<this->m; ++i){
            for(int j=0; j<other.n; ++j){
                double sum=0.0;
                for(int k=0; k<this->n; ++k){
                    sum+=this->elements[i][k]*other.elements[k][j];
                }
                newmatrix.elements[i][j]=sum;
            }
        }
    }
    else{
        throw invalid_argument("Matrix dimensions don't match for multiplication.");
    }
    return newmatrix;
}

Matrix Matrix::operator+(const Matrix &other) const{
    Matrix new_Matrix(this->m, this->n);
    if(this->m==other.m && this->n==other.n){
        for(int i=0; i<m; ++i){
            for(int j=0; j<n; ++j){
                new_Matrix.elements[i][j]=this->elements[i][j]+other.elements[i][j];
            }
        }
    }
    else{
        throw invalid_argument("Matrix dimensions don't match for plus.");
    }
    return new_Matrix;
}

void Matrix::add_elements(const Matrix &other){
    if(this->n == other.n){
        for(int i=0; i<other.m; ++i){
            elements.push_back(other.elements[i]);
        }
    }
    
}

void Matrix::set_elements(const Matrix &other,const int &index){
    if(this->n == other.n && other.m==0){
        this->elements[index]=other.elements[0];
    }
}

Matrix Matrix::operator-(const Matrix &other) const{
    Matrix new_Matrix(this->m, this->n);
    if(this->m==other.m && this->n==other.n){
        for(int i=0; i<m; ++i){
            for(int j=0; j<n; ++j){
                new_Matrix.elements[i][j]=this->elements[i][j]-other.elements[i][j];
            }
        }
    }
    else{
        throw invalid_argument("Matrix dimensions don't match for plus.");
    }
    return new_Matrix;
}

void Matrix::print() const{
    for(int i=0; i<this->m; ++i){
        for(int j=0; j<this->n; ++j){
            cout<<this->elements[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;
}