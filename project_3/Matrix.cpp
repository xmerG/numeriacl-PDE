#include"Matrix.h"

Matrix::Matrix(){
    m=0; 
    n=0;
    elements.clear();
    elements.shrink_to_fit();
}

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

Matrix Matrix::operator=(const Matrix &other){
    elements.clear();
    m=other.m;
    n=other.n;
    elements=other.elements;
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

Matrix Matrix::operator*(const double &value){
    for(int i=0; i<m; ++i){
        for(int j=0; j<n; ++j){
            elements[i][j]*=value;
        }
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

void Matrix::add_row_elements(const Matrix &other){
    if(this->n == other.n){
        for(int i=0; i<other.m; ++i){
            elements.push_back(other.elements[i]);
        }
        this->m+=other.m;
    }

    
}

void Matrix::add_col_elements(const Matrix& other){
    if(this->m=other.m){
        for(int i=0; i<other.n; ++i){
            for(int j=0; j<this->m; ++j){
                elements[j].push_back(other.elements[j][i]);
            }
        }
        this->n+=other.n;
    }
}

void Matrix::set_elements(const Matrix &other,const int &index,const int &i1, const int &i2, const int &j1, const int &j2){
    if(other.m==i2-i1+1 && other.n==j2-j1+1){
        for(int i=i1; i<=i2; ++i){
            for(int j=j1; j<=j2; ++j){
                elements[i][j]=other.elements[i-i1][j-j1];
            }
        }
    }
    else{
        throw out_of_range("Matrix indices out of range");
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

void Matrix::set_elements(const int &i, const int &j, const double &value){
    if(0<=i<this->m && 0<=j<this->n){
        elements[i][j]=value;
    }
    else{
        throw invalid_argument("Matrix dimensions don't match for plus.");
    }
}

//vector<Matrix> Matrix::LU(){

//}

//Matrix Matrix::operator/(const Matrix &other){
    
//}

double Matrix::l2_norm() const{
    double norm=-1.0;
    if(m==1 || n==1){
        norm=0.0;
        for(int i=0; i<m; ++i){
            for(int j=0; j<n; ++j){
                norm+=pow(elements[i][j], 2);
            }
        }
        return sqrt(norm);
    }
    else{
        cerr<<"not l2 norm for matrix"<<endl;
        return norm;
    }
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