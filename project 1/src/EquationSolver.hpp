#ifndef EQUATIONSOLVER
#define EQUATIONSOLVER

#include<iostream>
#include<vector>
#include<cmath>
#include <lapacke.h>
#include <nlohmann/json.hpp>
#include"Function.hpp"
#include"BoundaryCondition.hpp"
#include"Domain.hpp"
using namespace std;

template<typename D = Domain, typename BC= BoundaryCondition>
class EquationSolver{
private:
    vector<vector<double>> grids;  // denote the 2 dimension grids
    vector<double> values; 
    vector<double> BD_values; //
    int N;  //the number of grids in 1 dimension

    vector<vector<double>> coeffMatrix(){
        if(BC==BoundaryCondition::Dirichlet){
            vector<vector<double>> A(N-1,vector<double>(N-1,0));
            A[N-2][N-2]=-4.0;
            for(int i=0; i<N-2; ++i){
                A[i][i]=-4.0;
            //--------------------------------------------to be completed
            }
        }
    }
public:
    EquationSolver(){};

    //-\Delat=f, g denotes the boundary condition
    EquationSolver(const int &_N, const Function &f, const Function &g){  
        double h=1.0/N;
        if(BC==BoundaryCondition::Dirichlet){
            for(int j=1; j<N; ++j){
                for(int i=1; i<N; ++i){
                    values.push_back(f(i*h, j*h));
                    grids.push_back(vector{i*h,j*h});
                }
            }
        }
        else if(BC==BoundaryCondition::Neumann){
            for(int j=0; j<=N; ++j){
                for(int i=0; i<=N; ++i){
                    values.push_back(f(i*h, j*h));
                    grids.push_back(vector{i*h, j*h});
                }
            }
        }

    }
};


#endif