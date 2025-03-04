#ifndef EQUATIONSOLVER
#define EQUATIONSOLVER

#include<iostream>
#include<vector>
#include<cmath>
#include<fstream>
#include <lapacke.h>
#include <nlohmann/json.hpp>
#include"Function.hpp"
#include"BoundaryCondition.hpp"
#include"Domain.hpp"
using namespace std;

template<Domain D, BoundaryCondition BC>
class EquationSolver{
private:
    vector<vector<double>> grids;  // denote the 2 dimension grids
    vector<double> values; 
    int N=0;  //the number of grids in 1 dimension
    double h=0.0;

    vector<vector<double>> coeffMatrix(){
        vector<vector<double>> A((N-1)*(N-1),vector<double>((N-1)*(N-1),0.0));
        if(BC==BoundaryCondition::Dirichlet){
            for(int k=0; k<N-1; ++k){
                int m=k*(N-1);
                for(int i=0; i<N-1; ++i){
                    A[i+m][i+m]=4.0;
                    if(i!=N-2){
                        A[i+m][i+m+1]=-1.0;
                        A[i+m+1][i+m]=-1.0; 
                    }
                    if(k!=N-2){
                        A[m+i][m+i+N-1]=-1.0;
                        A[m+i+N-1][m+i]=-1.0;
                    }
                }
            }
        }
        else if(BC==BoundaryCondition::Neumann){
            //vector<vector<double>> A((N+1)*(N+1), vector<double>((N+1)*(N+1), 0.0));
            /*for(int k=0; k<N+1; ++k){
                int m=k*(N+1);
                if(k==0){
                    A[m][m]=2.0;
                    A[N+m][N+m]=2.0;
                    A[m][m+1]=-1.0;
                    A[m+N][m+N-1]=-1.0;
                    A[m][m+N+1]=-1.0;
                    A[m+N][m+2*N+1]=-1.0;
                    for(int i=1; i<N; ++i){
                        A[m+i][m+i]=4.0;
                        A[m+i][m+i-1]=-1.0;
                        A[m+i][m+i+1]=-1.0;
                        A[m+i][m+i+N+1]=-2.0;
                    }
                }
                else if(k==N){
                    A[m][m]=2.0;
                    A[N+m][N+m]=2.0;
                    A[m][m+1]=-1.0;
                    A[m+N][m+N-1]=-1.0;
                    A[m][m-N-1]=-1.0;
                    A[m+N][m-1]=-1.0;
                    for(int i=0; i<N; ++i){
                        A[m+i][m+i]=4.0;
                        A[m+i][m+i-1]=-1.0;
                        A[m+i][m+i+1]=-1.0;
                        A[m+i][m+i-N-1]=-2.0;
                    }
                }
                else{
                    for(int j=0; j<N+1; ++j){
                        A[m+j][m+j]=4.0;
                        A[m+j][m+j+N+1]=-1.0;
                        A[m+j][m+j-N-1]=-1.0;
                        if(j!=N && j!=0){
                            A[m+j][m+j+1]=-1.0;
                            A[m+j][m+j-1]=-1.0;
                        }
                        else if(j==0){
                            A[m][m+1]=-2.0;
                        }
                        else{
                            A[m+N][m+N-1]=-2.0;
                        }
                    }
                }
            }*/
            for(int k=0; k<N-1; ++k){
                int m=k*(N-1);
                if(k==0){
                    A[m][m]=1.0;
                    A[N-2+m][N-2+m]=4.0;
                    A[m+N-2][m+N-3]=-2.0;
                    A[m+N-2][m+2*N-3]=-2.0;
                    for(int i=1; i<N-2; ++i){
                        A[m+i][m+i]=8.0;
                        A[m+i][m+i-1]=-3.0;
                        A[m+i][m+i+1]=-3.0;
                        A[m+i][m+i+N-1]=-2.0;
                    }
                }
                else if(k==N-2){
                    A[m][m]=4.0;
                    A[N-2+m][N-2+m]=4.0;
                    A[m][m+1]=-2.0;
                    A[m+N-2][m+N-3]=-2.0;
                    A[m][m-N+1]=-2.0;
                    A[m+N-2][m-1]=-2.0;
                    for(int i=1
                        ; i<N-2; ++i){
                        A[m+i][m+i]=8.0;
                        A[m+i][m+i-1]=-3.0;
                        A[m+i][m+i+1]=-3.0;
                        A[m+i][m+i-N+1]=-2.0;
                    }
                }
                else{
                    for(int j=0; j<N-1; ++j){
                        if(j!=N-2 && j!=0){
                            A[m+j][m+j]=4.0;
                            A[m+j][m+j+1]=-1.0;
                            A[m+j][m+j-1]=-1.0;
                            A[m+j][m+j+N-1]=-1.0;
                            A[m+j][m+j-N+1]=-1.0;
                        }
                        else if(j==0){
                            A[m][m]=8.0;
                            A[m][m+1]=-2.0;
                            A[m][m+N-1]=-3.0;
                            A[m][m-N+1]=-3.0;
                        }
                        else{
                            A[m+N-2][m+N-2]=8.0;
                            A[m+N-2][m+N-3]=-2.0;
                            A[m+N-2][m+2*N-3]=-3.0;
                            A[m+N-2][m-1]=-3.0;
                        }
                    }
                }
            }
            // to be completed
            for(int i=0; i<(N-1)*(N-1); ++i){
                for(int j=0; j<(N-1)*(N-1); ++j){
                    cout<<A[i][j]<<" ";
                }
                cout<<endl;
            }
        }
        return A;
    }

    vector<double> convert(){
        vector<vector<double>> A=coeffMatrix();
        vector<double> a;
        for(int i=0; i<A.size(); ++i){
            for(int j=0; j<A.size(); ++j){
                a.push_back(A[j][i]);
            }
        }
        return a;
    }


    void getcolumn(const Function &g){
        int n=(N-1)*(N-1);
        if(BC==BoundaryCondition::Dirichlet){
            values[0]=values[0]+g(h, 0.0)+g(0.0, h);
            values[N-2]=values[N-2]+g(1.0,h)+g(1-h, 0.0);
            values[n-1]=values[n-1]+g(1.0, 1-h)+g(1-h,1.0);
            values[(N-2)*(N-1)]=values[(N-2)*(N-1)]+g(0.0,1-h)+g(h,1.0);
            for(int i=1; i<N-2; ++i){
                values[i]+=g((i+1)*h, 0.0);
                values[(N-1)*i]+=g(0.0, (i+1)*h);
                values[(N-1)*(i+1)-1]+=g(1.0, (i+1)*h);
                values[(N-1)*(N-2)+i]+=g((i+1)*h, 1.0);
            }
        }
        else if(BC==BoundaryCondition::Neumann){
            double m=2*h;
            values[0]=0.0;
            values[n-1]=3*values[n-1]+m*g(1-h, 1.0)+m*g(1.0,1-h);
            values[N-2]=3*values[N-2]+m*g(1-h,0.0)+m*g(1.0, h);
            values[(N-2)*(N-1)]=3*values[(N-2)*(N-1)]+m*g(h,1.0)+m*g(0.0,1-h);
            for(int i=0; i<N-1; ++i){
                int k=(N-1)*i;
                    if(i==0){
                        for(int j=1; j<N-2; ++j){
                            values[j]=3*values[j]+m*g((j+1)*h, 0.0);
                        }
                    }
                    else if(i==N-2){
                        for(int j=1; j<N-2; ++j){
                            values[k+j]=3*values[k+j]+m*g((j+1)*h, 1.0);
                        }
                    }
                    else{
                        values[k]=3*values[k]+m*g(0.0, (i+1)*h);
                        values[k+N-2]=3*values[k+N-2]+m*g(1.0, (i+1)*h);
                    }
            }
            //-----------------------------------------------------------------
            
            for(int i=0; i<n;++i){
                cout<<values[i]<<"  ";
            }
            cout<<endl;
        }
    }
    void solve(const Function &g){
        vector<double> matrix=convert();
        int n=(N-1)*(N-1);
        getcolumn(g);
        vector<int> ipiv(n);
        int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, n, 1, matrix.data(), n, ipiv.data(), values.data(), n);

    }
public:
    EquationSolver(){};

    //-\Delat=f, g denotes the boundary condition
    EquationSolver(const int &_N, const Function &f){  
        N=_N;
        h=1.0/N;
        for(int j=1; j<N; ++j){
            for(int i=1; i<N; ++i){
                values.push_back(f(i*h, j*h)*h*h);
                grids.push_back(vector{i*h,j*h});
            }
        }
    }

    vector<vector<double>> getgrids(){
        return grids;
    }

    vector<double> getvalues(){
        return values;
    }

    vector<double> errors(const Function &f){
        int n=grids.size();
        vector<double> e(n, 0.0);
        for(int i=0; i<n; ++i){
            double x=grids[i][0];
            double y=grids[i][1];
            e[i]=abs(values[i]-f(x,y));
        }
        return e;

    }

    void norm_error(const Function &f){
        vector<double> error=this->errors(f);
        double l1_norm=0.0;
        double l2_norm=0.0;
        double infinity_norm=0.0;
        for(int i=0; i<error.size(); ++i){
            l1_norm+=h*error[i];
            l2_norm+=h*pow(error[i], 2);
            if(error[i]>infinity_norm){
                infinity_norm=error[i];
            }
        }
        cout<<"------------------------------------------- errors -----------------------------------------"<<endl;
        cout<<"l_1 norm is "<<l1_norm<<endl;
        cout<<"l_2 norm is "<<sqrt(l2_norm)<<endl;
        cout<<"l_infty norm is"<<infinity_norm<<endl;

    }

    void solveEquation(const Function &g){
        solve(g);
    }

    void print(const string &filename, const Function &f){
        nlohmann::json j;
        j["boundary_condition"] = BC; 
        j["grids"] = grids;
        j["values_on_grids"] = values;
        vector<double> errors = this->errors(f);
        j["errors"] = errors;  
        std::ifstream file_check(filename); 
        bool is_empty = file_check.peek() == std::ifstream::traits_type::eof(); 
        file_check.close();  
        nlohmann::json jsonDataArray;
        if (!is_empty) {
            std::ifstream inFile(filename);  
            inFile >> jsonDataArray;  
            inFile.close();  
        }
        jsonDataArray.push_back(j);
        std::ofstream outFile(filename, std::ios::out | std::ios::trunc);  // 打开文件，清空内容
        if (outFile.is_open()) {
            // 将修改后的 JSON 数组写回文件，并格式化输出
            outFile << jsonDataArray.dump(4);  // 4 个空格缩进
            outFile.close();
        } 
        else {
            cerr << "Error opening file " << filename << endl;
        }
    }
    
};


#endif