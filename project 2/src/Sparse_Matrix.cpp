#include"Sparse_Matrix.h"
Sparse_Matrix::Sparse_Matrix():n(0){}

Sparse_Matrix::Sparse_Matrix(const int &_n):n(_n){}

Sparse_Matrix::Sparse_Matrix(const int &_n, const vector<label> &e):n(_n), elements(e){}

void Sparse_Matrix::setValues(const int &i, const int &j, const double &value){
    if(i>=0 && i<n && j>=0 && j<n){
        elements[i][j]=value;
    }
    else{
        cerr<<"invalid label!"<<endl;
    }
}

double Sparse_Matrix::operator()(const int &i, const int &j) {
    if(i>=0 && i<n && j>=0 && j<n){
        return elements[i][j];
    }
    else{return 0.0;}
}

//void Sparse_Matrix::transform(){
    //for(map<label, double>::iterator i=A.begin(); i!=A.end(); ++i){
        
    //}
//}

Sparse_Matrix Sparse_Matrix::operator+(const Sparse_Matrix &B){
    if(n!=B.n){
        cerr<<"invalid since the matrix dimension must agree for addision"<<endl;
        return Sparse_Matrix();
    } 
    else{
        vector<label> new_elements(n);
        for(int i=0; i<n; ++i){
            label lA=elements[i];
            label l=B.elements[i];
            new_elements[i]=lA;
            for(map<int, double>::iterator itr=l.begin(); itr!=l.end(); ++itr){
                if(lA.count(l[itr->first])){
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

Sparse_Matrix Sparse_Matrix::operator*(const double &a){
    vector<label> new_elements(n);
    for(int i=0; i<n; ++i){
        label l=elements[i];
        new_elements[i]=elements[i];
        for(map<int,double>::iterator itr=l.begin(); itr!=l.end(); ++itr){
            new_elements[i][itr->first]=a*itr->second;
        }
    }
    return Sparse_Matrix(n, new_elements);
}

Vector Sparse_Matrix::operator*(const Vector &v){
    vector<double> new_elements(n);
    vector<double> V=v.elements;
    for(int i=0; i<n; ++i){
        label l=elements[i];
        double current=0.0;
        for(map<int, double>::iterator itr=l.begin(); itr!=l.end(); ++itr){
            current+=itr->second*V[itr->first];
        }
        new_elements[i]=current;
    }
    return Vector(n, new_elements);
}