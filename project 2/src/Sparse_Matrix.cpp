#include"Sparse_Matrix.h"
Sparse_Matrix::Sparse_Matrix():n(0){}

Sparse_Matrix::Sparse_Matrix(const int &_n):n(_n){
    elements.resize(n);
}

Sparse_Matrix::Sparse_Matrix(const int &_n, const vector<label> &e):n(_n), elements(e){}

/*Sparse_Matrix::Sparse_Matrix(Sparse_Matrix&& other) noexcept: n(other.n), elements(move(other.elements)){
    other.n=0;
}*/

/*Sparse_Matrix& Sparse_Matrix::operator=(Sparse_Matrix&& other) noexcept{
    if (this != &other) {
        n = other.n;
        elements = move(other.elements);
        other.n = 0;
    }
    return *this;
}*/

void Sparse_Matrix::setValues(const int &i, const int &j, const double &value){
    if(i>=0 && i<n && j>=0 && j<n){
        label l=elements[i];
        l[j]=value;
        elements[i]=l;
    }
}

double Sparse_Matrix::operator()(const int &i, const int &j) const{
    if(i>=0 && i<n && j>=0 && j<n){
        map<int, double>::const_iterator itr=elements[i].find(j);
        if(itr!=elements[i].end()){
            return itr->second;
        }
        else{return 0.0;}
    }
    else{return 0.0;}
}

//void Sparse_Matrix::transform(){
    //for(map<label, double>::iterator i=A.begin(); i!=A.end(); ++i){
        
    //}
//}

Sparse_Matrix Sparse_Matrix::operator+(const Sparse_Matrix &B) const{
    if(n!=B.n){
        cerr<<"invalid since the matrix dimension must agree for addision"<<endl;
        return Sparse_Matrix();
    } 
    else{
        vector<label> new_elements(n);
        for(int i=0; i<n; ++i){
            new_elements[i]=elements[i];
            label l=B.elements[i];
            for(map<int, double>::const_iterator itr=l.begin(); itr!=l.end(); ++itr){
                if(new_elements[i].count(l[itr->first])){
                    new_elements[i][itr->first]+=itr->second;
                }
                else{
                    new_elements[i][itr->first]=itr->second;
                }
            }
        }
        return Sparse_Matrix(n, new_elements);
    }
}

Sparse_Matrix Sparse_Matrix::operator*(const double &a) const{
    vector<label> new_elements(n);
    for(int i=0; i<n; ++i){
        for(map<int,double>::const_iterator itr=elements[i].begin(); itr!=elements[i].end(); ++itr){
            new_elements[i][itr->first]=a*itr->second;
        }
    }
    return Sparse_Matrix(n, new_elements);
}

Vector Sparse_Matrix::operator*(const Vector &v) const{
    vector<double> new_elements(n);
    vector<double> V=v.elements;
    for(int i=0; i<n; ++i){
        double current=0.0;
        for(map<int, double>::const_iterator itr=elements[i].begin(); itr!=elements[i].end(); ++itr){
            current+=itr->second*V[itr->first];
        }
        new_elements[i]=current;
    }
    return Vector(n, new_elements);
}

vector<double> Sparse_Matrix::convert_to_vector() const{
    int size=n*n;
    vector<double> v(size, 0.0);
    for(int i=0; i<n; ++i){
        for(map<int, double>::const_iterator itr=elements[i].begin(); itr!=elements[i].end(); ++itr){
            v[n*i+itr->first]=itr->second;
        }
    }
    return v;
}

int Sparse_Matrix::getdim() const{return n;}

void Sparse_Matrix::print() {
    for(int i=0; i<n; ++i){
        cout<<"[";
        for(int j=0; j<n; ++j){
            if(elements[i].count(j)){
                cout<<elements[i][j]<<", ";
            }
            else{
                cout<<"0"<<", ";
            }
        }
        cout<<"],"<<endl;
    }
}

void Sparse_Matrix::Gauss_Seidel(Vector &initial, const Vector &b){
    for(int i=0; i<n; ++i){
        label l=elements[i];
        double value=0.0;
        double a=0.0;
        for(map<int, double>::iterator itr=l.begin(); itr!=l.end(); ++itr){
            if(itr->first!=i){
                value-=itr->second*initial(itr->first);
            }
            else{
                a=itr->second;
            }
        }
        value+=b(i);
        value/=a;
        initial.set_Value(i, value);
    }
}
/*
void Sparse_Matrix::solve(Vector b){
    vector<double> B=b.getelements();
    vector<double> matrix = this->convert_to_vector();  
    vector<int> ipiv(n); 
    int info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, 1, matrix.data(), n, ipiv.data(), B.data(), n);
    b=Vector(n, B);
}*/
