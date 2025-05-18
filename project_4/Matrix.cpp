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

int Matrix::get_col() const{return n;}

int Matrix::get_row() const{return m;}

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

Matrix Matrix::operator()(const int &i1, const int &i2, const int &j1, const int &j2){
    vector<vector<double>> new_element(i2-i1+1, vector<double>(j2-j1+1, 0.0));
    if(0<=i1 && i1<=i2 && i2<=this->m && 0<=j1 && j1<=j2 && j2<=this->n){
        for(int i=0; i<=i2-i1; ++i){
            for(int j=0; j<=j2-j1; ++j){
                new_element[i][j]=this->elements[i+i1][j+j1];
            }
        }
    }
    return move(Matrix(i2-i1+1, j2-j1+1, new_element));
}

Matrix Matrix::operator*(const double &value){
    vector<vector<double>> new_elements(m,vector<double>(n, 0.0));
    for(int i=0; i<m; ++i){
        for(int j=0; j<n; ++j){
            new_elements[i][j]=value*elements[i][j];
        }
    }
    return move(Matrix(m, n, new_elements));
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
    return move(newmatrix);
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
    return move(new_Matrix);
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

void Matrix::set_row_elements(const Matrix &other, const int &index){
    if(this->n == other.n){
        this->elements[index]=other.elements[0];
    }
}

void Matrix::set_col_elements(const Matrix &other, const int &index){
    if(this->m == other.m){
        for(int i=0; i<m; ++i){
            this->elements[i][index]=other.elements[i][0];
        }
    }
}

void Matrix::set_elements(const Matrix &other,const int &i1, const int &i2, const int &j1, const int &j2){
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

vector<Matrix> Matrix::LU() const{
    vector<vector<double>> L_elements(m, vector<double>(m, 0.0));
    vector<vector<double>> U_elements(n, vector<double>(n, 0.0));
    vector<vector<double>> p_elements(n, vector<double>(n, 0.0));
    for(int i=0; i<n; ++i){
        p_elements[i][i]=1.0;
    }
    vector<vector<double>> temp=elements;
    if(m==n){
        for(int i=0; i<n; ++i){
            int index=i;
            double value=abs(temp[i][i]);
            for(int j=i; j<n; ++j){
                if(value<abs(temp[j][i])){
                    value=abs(temp[i][j]);
                    index=j;
                }
            }
            vector<double> old=temp[index];
            temp[index]=temp[i];
            temp[i]=old;
            old=p_elements[index];
            p_elements[index]=p_elements[i];
            p_elements[i]=move(old);
            for(int j=i; j<n; ++j){
                L_elements[j][i]=temp[j][i]/temp[i][i];
                U_elements[i][j]=temp[i][j];
                if(j!=i){
                    for(int k=i+1; k<n; ++k){
                        temp[j][k]-=L_elements[j][i]*U_elements[i][k];
                    }
                }
            }
        }
    }
    Matrix L(m, m, L_elements);
    Matrix U(n, n ,U_elements);
    Matrix P(n, n, p_elements);
    vector<Matrix> LUP{move(L), move(U), move(P)};
    return LUP;
}

Matrix Matrix::operator/(Matrix &B){
    vector<Matrix> LUP=this->LU();
    B=move(LUP[2]*B);
    if(m==n && this->m==B.m){
        for(int k=0; k<B.m; ++k){
            for(int j=0; j<m-1; ++j){
                B.set_elements(j, k, B(j, k)/LUP[0](j, j));
                for(int i=j+1; i<m; ++i){
                    B.set_elements(i, k, B(i, k)-B(i, k)*LUP[0](i, j));
                }
            }
            B.set_elements(m-1, k, B(m-1, k)/LUP[0](m-1, m-1));
            for(int j=m-1; j>0; --j){
                B.set_elements(j, k, B(j, k)/LUP[1](j, j));
                for(int i=1; i<j-1; ++i){
                    B.set_elements(i, k, B(i, k)-B(i, k)*LUP[1](i, j));
                }
            }
            B.set_elements(0,k, B(0,k)/LUP[1](0,0));
        }
    }
    return move(B);

}

Matrix Matrix::sparse_cholesky() const{
    vector<vector<double>> new_elements(m, vector<double>(n, 0.0));
    if(this->m==this->n){
        new_elements[0][0]=sqrt(elements[0][0]);
        new_elements[1][0]=elements[1][0]/new_elements[0][0];
        for(int i=1; i<this->m-1; ++i){
            new_elements[i][i]=sqrt(elements[i][i]-pow(new_elements[i][i-1], 2));
            new_elements[i+1][i]=elements[i+1][i]/new_elements[i][i];
        }
        new_elements[n-1][n-1]=sqrt(elements[n-1][n-1]-pow(new_elements[n-1][n-2], 2));
    }
    return move(Matrix(m, n, new_elements));
}

Matrix Matrix::solve(const Matrix &b) const{
    Matrix l=move(this->sparse_cholesky());
    vector<vector<double>> sol(n, vector<double>(1, 0.0));
    sol[n-1][0]=b.elements[n-1][0]/l.elements[n-1][n-1];
    for(int i=m-2; i>=0; --i){
        sol[i][0]=(b.elements[i][0]-l.elements[i+1][i]*sol[i+1][0])/l.elements[i][i];
    }
    sol[0][0]/=l.elements[0][0];
    for(int i=1; i<m; ++i){
        sol[i][0]=(sol[i][0]-sol[i+1][0]*l.elements[i][i-1])/l.elements[i][i];
    }
    return move(Matrix(n, 1, sol));
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