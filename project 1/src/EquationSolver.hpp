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
#include"Circle.hpp"
#include"Domain.hpp"
using namespace std;

template<Domain D, BoundaryCondition BC>
class EquationSolver{
private:
    vector<vector<double>> grids;  // denote the 2 dimension grids
    vector<double> values; 
    int N=0;  //the number of grids in 1 dimension
    double h=0.0;
    Circle *c=nullptr;
    vector<bool> incircle;
    vector<double> laplacian;
    vector<vector<double>> A;

    void precondition(){
        A=vector<vector<double>>((N-1)*(N-1),vector<double>((N-1)*(N-1),0.0));
        for(int k=1; k<N-2; ++k){
            int m=k*(N-1);
            for(int i=1; i<N-2; ++i){
                int index=m+i;
                A[index][index]=4.0;
                A[index][index-1]=-1.0;
                A[index][index+1]=-1.0;
                A[index][index+N-1]=-1.0;
                A[index][index-N+1]=-1.0;
            }
        }
    }


    void coeffMatrix(const vector<int> & mixed){
        A=vector<vector<double>>((N-1)*(N-1),vector<double>((N-1)*(N-1),0.0));
        if(BC==BoundaryCondition::Neumann){
            /*this->precondition();
            for(int i=1; i<N-1; ++i){
                A[i][i]=8.0;
                A[i][i+N-1]=-2.0;
                A[i][i-1]=-3.0;
                A[i][i+1]=-3.0;
                int index=i*(N-1);
                A[index][index]=8.0;
                A[index][index+1]=-2.0;
                A[index][index+N-1]=-3.0;
                A[index][index-N+1]=-3.0;
                index+=N-2;
                A[index][index]=8.0;
                A[index][index-1]=-2.0;
                A[index][index+N-1]=-3.0;
                A[index][index-N+1]=-3.0;
                index=N*N-3*N+1+i;
                A[index][index]=8.0;
                A[index][index+1]=-3.0;
                A[index][index-1]=-3.0;
                A[index][index-N+1]=-2.0;
            }
            A[0][0]=1.0;
            A[N-2][N-2]=4.0;
            A[N-2][N-3]=-2.0;
            A[N-2][2*N-3]=-2.0;
            A[(N-1)*(N-1)-1][(N-1)*(N-1)-1]=4.0;
            A[(N-1)*(N-1)-1][N*N-2*N-1]=-2.0;
            A[(N-1)*(N-1)-1][N*N-3*N+1]=-2.0;
            A[(N-1)*(N-2)][(N-1)*(N-2)]=4.0;
            A[(N-1)*(N-2)][(N-1)*(N-2)+1]=-2.0;
            A[(N-1)*(N-2)][(N-1)*(N-3)]=-2.0;
            for(int i=0;i<(N-1)*(N-1); ++i){
                for(int j=0; j<(N-1)*(N-1); ++j){
                    cout<<A[i][j]<<" ";
                }
                cout<<endl;
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
        }
        else{
            this->precondition();
            for(int i=0; i<N-1; ++i){
                bool p=true;
                if(i==0 || i==N-2){
                    p=false;
                }
                int index=i*(N-1);
                if(mixed[0]==0){
                    A[i][i+N-1]=-1.0;
                    if(p){
                        A[i][i]=4.0;
                        A[i][i-1]=-1.0;
                        A[i][i+1]=-1.0;
                    }
                    else{
                        A[i][i]+=2.0;
                    }
                }
                else{
                    if(p){
                        A[i][i]=8.0;
                        A[i][i+1]=-3.0;
                        A[i][i-1]=-3.0;
                        A[i][i+N-1]=-2.0;
                    }
                    else{
                        A[i][i]+=2.0/3.0;
                        A[i][i+N-1]=-2.0/3.0;
                    }
                }

                if(mixed[1]==0){
                    A[index][index+1]=-1.0;
                    if(p){
                        A[index][index]=4.0;
                        A[index][index+N-1]=-1.0;
                        A[index][index-N+1]=-1.0;
                    }
                    else{
                        A[index][index]+=2.0;
                        
                    }
                }
                else{
                    if(p){
                        A[index][index]=8.0;
                        A[index][index+1]=-2.0;
                        A[index][index+N-1]=-3.0;
                        A[index][index-N+1]=-3.0;
                    }
                    else{
                        A[index][index]+=2.0/3.0;
                        A[index][index+1]=-2.0/3.0;
                    }
                }
                index+=N-2;
                if(mixed[2]==0){
                    A[index][index-1]=-1.0;
                    if(p){
                        A[index][index]=4.0;
                        A[index][index+N-1]=-1.0;
                        A[index][index-N+1]=-1.0;
                    }
                    else{
                        A[index][index]+=2.0;
                    }
                }
                else{
                    if(p){
                        A[index][index]=8.0;
                        A[index][index-1]=-2.0;
                        A[index][index+N-1]=-3.0;
                        A[index][index-N+1]=-3.0;
                    }
                    else{
                        A[index][index]+=2.0/3.0;
                        A[index][index-1]=-2.0/3.0;
                    }
                }
                index=(N-1)*(N-2)+i;
                if(mixed[3]==0){
                    A[index][index-N+1]=-1.0;
                    if(p){
                        A[index][index]=4.0;
                        A[index][index+1]=-1.0;
                        A[index][index-1]=-1.0;
                    }
                    else{
                        A[index][index]+=2.0;
                    }
                }
                else{
                    if(p){
                        A[index][index]=8.0;
                        A[index][index-N+1]=-2.0;
                        A[index][index-1]=-3.0;
                        A[index][index+1]=-3.0;
                    }
                    else{
                        A[index][index]+=2.0/3.0;
                        A[index][index-N+1]=-2.0/3.0;
                    }
                }
            }
            /*if(BC==BoundaryCondition::Neumann){
                A[0]=vector<double>((N-1)*(N-1), 0.0);
                A[0][0]=1.0;
            }*/
        }
    }

    void Diri(const Function &g){
        for(int k=0; k<N-1; ++k){
            int m=k*(N-1);
            for(int i=0; i<N-1; ++i){
                int index=i+m;
                double current_x=grids[index][0];
                double current_y=grids[index][1];
                if(incircle[index]==true){
                    A[index]=vector<double>((N-1)*(N-1), 0.0);
                    A[index][index]=1.0;
                }
                else{
                    double dx=c->x_distance_to_circle(current_x, current_y);
                    double dy=c->y_distance_to_circle(current_x, current_y);
                    double alpha=1.0;
                    double theta=1.0;
                    if(abs(dx)<=h && abs(dx)>0){
                        int direction=1;
                        if(dx<0){
                            direction=-1;
                        }
                        theta=abs(dx)/h;
                        values[index]+=2.0*g((i+1)*h+dx, (k+1)*h)/(theta*(1+theta));
                        A[index][index+direction]=0.0;
                        if(i!=0 && i!=N-2){//左右两边都不是正方形边界
                            A[index][index-direction]=-2.0/(1+theta);
                        }
                        else if(i==0){ //左边是正方形边界
                            values[index]+=2.0*g(0.0, (k+1)*h)/(1+theta);
                        }
                        else{
                            values[index]+=2.0*g(1.0,(k+1)*h)/(1+theta);
                        }
                    }
                    else{
                        if(i==0){
                            values[index]+=g(0.0, (k+1)*h);
                        }
                        else if(i==N-2){
                            values[index]+=g(1.0, (k+1)*h);
                        }
                    }
                    if(abs(dy)<=h && abs(dy)>0){
                        int direction=1;
                        if(dy<0){
                            direction=-1;
                        }
                        alpha=abs(dy)/h;
                        values[index]+=2.0*g((i+1)*h, (k+1)*h+dy)/(alpha*(1+alpha));
                        A[index][index+(N-1)*direction]=0.0;
                        if(k!=0 && k!=N-2){
                            A[index][index-(N-1)*direction]=-2.0/(1+alpha);
                        }
                        else if(k==0){//下面是正方形边界
                            values[index]+=2.0*g((i+1)*h, 0.0)/(1+alpha);
                        }
                        else{
                            values[index]+=2.0*g((i+1)*h, 1.0)/(1+alpha);
                        }
                    }
                    else{
                        if(k==0){
                            values[index]+=g((i+1)*h, 0.0);
                        }
                        else if(k==N-2){
                            values[index]+=g((i+1)*h, 1.0);
                        }
                    }
                    A[index][index]=2.0/alpha+2.0/theta;
                }
            }
        }
    }

    void Neum(const Function &g, const vector<int> &mixed){
        this->getcolumn(g,mixed);
        for(int k=0; k<N-1; ++k){
            int m=k*(N-1);
            for(int i=0; i<N-1; ++i){
                int index=i+m;
                double current_x=grids[index][0];
                double current_y=grids[index][1];
                if(incircle[index]==true){
                    A[index]=vector<double>((N-1)*(N-1), 0.0);
                    A[index][index]=1.0;
                    values[index]=0.0;
                }
                else{
                    double dx=c->x_distance_to_circle(current_x, current_y);
                    double dy=c->y_distance_to_circle(current_x, current_y);
                    int xdirection=1;
                    int ydirection=1;
                    double r=c->get_radius();
                    bool p=(i==0 || i==N-2 || k==0 || k==N-2);
                    if(abs(dx)<h && abs(dx)>0){                            
                        if(dx<0){
                            xdirection=-1;
                        }
                        A[index][index+xdirection]=0.0;
                        
                        if(dy==0.0){
                            if(p){
                                A[index][index]-=2.0;;
                                values[index]+=2.0*h*g(dx+current_x,current_y);
                            }
                            else{
                                A[index][index]-=1.0;;
                                values[index]+=h*g(current_x+dx,current_y);
                            }
                        }
                        else{
                            if(dy<0){ydirection=-1;}
                            double temp=(current_y-c->getY())/(abs(current_x-c->getX())-h);
                            double Ty=current_y+h*temp;
                            double d=c->distance(current_x, Ty);
                            double Ex=current_x+(c->getX()-current_x)*(d-r)/d;
                            double Ey=Ty+(c->getY()-Ty)*(d-r)/d;
                            if(p){
                                A[index][index]+=-2.0+2.0*abs(temp);
                                A[index][index-(N-1)*ydirection]+=-2.0*abs(temp);
                                values[index]+=2.0*h*d*g(Ex, Ey)/(abs(current_x-c->getX()));
                            }
                            else{
                                A[index][index]+=-1.0+abs(temp);
                                A[index][index-(N-1)*ydirection]+=-abs(temp);
                                values[index]+=h*d*g(Ex, Ey)/(abs(current_x-c->getX()));
                            }
                        }                    
                    }
                    if(abs(dy)<h && abs(dy)>0){
                        if(dy<0){
                            ydirection=-1;
                        }
                        A[index][index+(N-1)*ydirection]=0.0;
                        if(dx==0.0){
                            if(p){
                                A[index][index]-=2.0;;
                                double temp=c->getY()-ydirection*r;
                                values[index]+=2.0*h*g(current_x, temp);
                            }
                            else{
                                A[index][index]-=1.0;;
                                double temp=c->getY()-ydirection*r;
                                values[index]+=h*g(current_x, temp);
                            }
                        }
                        else{
                            if(dx<0){
                                xdirection=-1;
                            }
                            double temp=(current_x-c->getX())/(abs(current_y-c->getY())-h);
                            double Tx=current_x+h*temp;
                            double d=c->distance(Tx, current_y);
                            double Ex=Tx+(c->getX()-Tx)*(d-r)/d;
                            double Ey=current_y+(c->getY()-current_y)*(d-r)/d;
                            if(p){
                                A[index][index]+=-2.0+2.0*abs(temp);
                                A[index][index-xdirection]+=-2.0*abs(temp);
                                values[index]+=2.0*d*h*g(Ex, Ey)/abs(current_y-c->getY());
                            }
                            else{
                                A[index][index]+=-1.0+abs(temp);
                                A[index][index-xdirection]+=-abs(temp);
                                values[index]+=d*h*g(Ex, Ey)/abs(current_y-c->getY());
                            }
                        }
                    }
                }
            }
        }
    }
    void coeffMatrix(const Function &g, const vector<int> &mixed){
        if(BC==BoundaryCondition::Dirichlet){
            this->coeffMatrix(mixed);
            this->Diri(g);
        }
        else if(BC==BoundaryCondition::Neumann){
            this->coeffMatrix(mixed);
            this->Neum(g, mixed);
        }
        else{
            this->coeffMatrix(mixed);
            if(mixed[4]==1){
                this->Neum(g,mixed);
            }
            else{
                this->Diri(g);
            }
        }

    }


    /*vector<vector<double>> coeff_Matrix(const Function &g, const vector<double> &Diri){
        vector<vector<double>> A((N+1)*(N+1),vector<double>((N+1)*(N+1),0.0));
        for(int k=0; k<N+1; ++k){
            int m=k*(N+1);
            for(int i=0; i<=N; ++i){
                int index=m+i;
                if(incircle[index]){
                    A[index][index]=1.0;
                    values[index]=0.0;
                }
                else{
                    if(k==0){
                        if(i!=0 && i!=N){
                            A[index][index+N+1]=-2.0;
                            A[index][index-1]=-1.0;
                            A[index][index+1]=-1.0;
                            values[index]+=2*h*g(i*h, 0.0);
                            A[index][index]=4.0;
                        }
                        else{
                            A[index][index]=1.0;
                        }

                    }
                    else if(k==N){
                        if(i!=0 &&i!=N){
                            A[index][index-1]=-1.0;
                            A[index][index+1]=-1.0;
                            A[index][index-N-1]=-2.0;
                            values[index]+=2*h*g(i*h, 1.0);
                            A[index][index]=4.0;
                        }
                        else{
                            A[index][index]=1.0;
                        }
                    }
                    else{
                        if(i==0){
                            A[index][index+1]=-2.0;
                            A[index][index]=4.0;
                            A[index][index+N+1]=-1.0;
                            A[index][index-N-1]=-1.0;
                            values[index]+=2*h*g(0.0,k*h);
                        }
                        else if(i==N){
                            A[index][index-1]=-2.0;
                            A[index][index]=4.0;
                            A[index][index+N+1]=-1.0;
                            A[index][index-N-1]=-1.0;
                            values[index]+=2*h*g(1.0, k*h);
                        }
                        else{
                            double current_x=grids[index][0];
                            double current_y=grids[index][1];
                            double dx=c->x_distance_to_circle(current_x, current_y);
                            double dy=c->y_distance_to_circle(current_x, current_y);
                            double alpha=1.0;
                            double theta=1.0;
                            if(abs(dx)<h){
                                theta=abs(dx)/h;
                                A[index][index]+=-2.0/(2*theta-1);
                                int direction=1;
                                if(dx<0){
                                    direction=-1;
                                }
                                A[index][index-direction]=2.0/(2*theta-1);
                                values[index]+=2.0*h*g(i*h+dx, k*h)*c->angle_x_direction(i*h+dx)/(2*theta-1);
                            }
                            else{
                                A[index][index]+=2.0;
                                A[index][index+1]=-1.0;
                                A[index][index-1]=-1.0;
                            }
                            if(abs(dy)<h){
                                alpha=abs(dy)/h;
                                int direction=1;
                                A[index][index]+=-2.0/(2*alpha-1);
                                if(dy<0){
                                    direction=-1;
                                }
                                A[index][index-direction*(N+1)]=2.0/(2*alpha-1);
                                values[index]+=2.0*h*g(i*h, k*h+dy)*c->angle_y_direction(k*h+dy)/(2*alpha-1);
                            }
                            else{
                                A[index][index]+=2.0;
                                A[index][index+N+1]=-1.0;
                                A[index][index-N-1]=-1.0;
                            }
                        }
                    }
                }

            }
        }
        values[0]=Diri[0];
        values[N]=Diri[1];
        values[(N+1)*N]=Diri[2];
        values[(N+1)*(N+1)-1]=Diri[3];
        for(int i=0; i<(N+1)*(N+1); ++i){
            cout<<values[i]<<" ";
        }
        cout<<endl;
        return A;
    }*/

    vector<double> convert(const Function &g, const vector<int> &mixed){
        if(D==Domain::regular){
            coeffMatrix(mixed);
        }
        else{
            coeffMatrix(g, mixed);
        }
        vector<double> a;
        for(int i=0; i<A.size(); ++i){
            for(int j=0; j<A.size(); ++j){
                a.push_back(A[j][i]);
            }
        }
        return a;
        A.clear();
        A.shrink_to_fit();
    }


    void getcolumn(const Function &g, const vector<int> &mixed){
        int n=(N-1)*(N-1);
        if(BC==BoundaryCondition::Neumann){
            double m=2*h;
            values[0]=0.0;
            values[n-1]=3*values[n-1]+m*g(1-h, 1.0)+m*g(1.0,1-h);
            values[N-2]=3*values[N-2]+m*g(1-h,0.0)+m*g(1.0, h);  
            values[(N-2)*(N-1)]=3*values[(N-2)*(N-1)]+m*g(h,1.0)+m*g(0.0,1-h);
            for(int i=1; i<N-2; ++i){
                values[i]=3*values[i]+m*g((i+1)*h, 0.0);
                values[i*(N-1)]=3*values[i*(N-1)]+m*g(0.0,(i+1)*h);
                values[(i+1)*(N-1)-1]=3*values[(i+1)*(N-1)-1]+m*g(1.0,(i+1)*h);
                values[(N-2)*(N-1)+i]=3*values[(N-2)*(N-1)+i]+m*g((i+1)*h, 1.0);
            }
            /*for(int i=0; i<N-1; ++i){
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
            }*/
        }
        else{
            if(mixed[0]==0){
                for(int i=0; i<N-1; ++i){
                    values[i]+=g((i+1)*h, 0.0);
                }
            }
            else{
                double m=2.0*h;
                for(int i=1; i<N-2; ++i){
                    values[i]=3*values[i]+m*g((i+1)*h,0.0);
                }
                values[0]+=m*g(h,0.0)/3.0;
                values[N-2]+=m*g(1-h,0.0)/3.0;
            }

            if(mixed[1]==0){
                for(int i=0; i<N-1; ++i){
                    values[(N-1)*i]+=g(0.0, (i+1)*h);
                }
            }
            else{
                double m=2.0*h;
                for(int i=1; i<N-2; ++i){
                    values[i*(N-1)]=3*values[i*(N-1)]+m*g(0.0, (i+1)*h);
                }
                values[0]+=m*g(0.0, h)/3.0;
                values[(N-2)*(N-1)]+=m*g(0.0, 1.0-h)/3.0;
            }

            if(mixed[2]==0){
                for(int i=0; i<N-1; ++i){
                    values[(i+1)*(N-1)-1]+=g(1.0, (i+1)*h);
                }
            }
            else{
                double m=2.0*h;
                for(int i=1; i<N-2; ++i){
                    values[(i+1)*(N-1)-1]=3*values[(i+1)*(N-1)-1]+m*g(1.0, (i+1)*h);
                }
                values[N-2]+=g(1.0,h)*m/3.0;
                values[N*N-2*N]+=g(1.0,1.0-h)*m/3.0;
            }

            if(mixed[3]==0){
                for(int i=0; i<N-1; ++i){
                    values[(N-1)*(N-2)+i]+=g((i+1)*h, 1.0);
                }
            }
            else{
                double m=2.0*h;
                for(int i=1; i<N-2; ++i){
                    values[(N-1)*(N-2)+i]=3*values[(N-1)*(N-2)+i]+m*g((i+1)*h, 1.0);
                }
                values[(N-1)*(N-2)]+=m*g(h,1.0)/3.0;
                values[N*N-2*N]+=m*g(1.0-h,1.0)/3.0;
            }
            /*if(BC==BoundaryCondition::Neumann){
                values[0]=0.0;
            }*/
            cout<<endl;
        }
    }
    void solve(const Function &g, const vector<int> &mixed){
        vector<double> matrix=convert(g,mixed);
        int n=(N-1)*(N-1);
        if(D==Domain::regular){
            getcolumn(g,mixed);
        }
        vector<int> ipiv(n);
        int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, n, 1, matrix.data(), n, ipiv.data(), values.data(), n);

    }

    void adjust(const double &initial){
        if(D==Domain::irregular){
            vector<vector<double>> newgrid;
            vector<double> newvalues;
            for(int i=0; i<(N-1)*(N-1); ++i){
                if(!incircle[i]){
                    newgrid.push_back(grids[i]);
                    newvalues.push_back(values[i]);
                }
            }
            grids=newgrid;
            values=newvalues;
        }
        if(BC==BoundaryCondition::Neumann){
            for(int i=0; i<values.size(); ++i){
                values[i]+=initial;
            }
        }
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

    vector<double> realValues(const Function &f){
        vector<double> realvalues;
        for(int i=0; i<grids.size(); ++i){
            double x=grids[i][0];
            double y=grids[i][1];
            realvalues.push_back(f(x,y));
        }
        return realvalues;
    }
public:
    EquationSolver(){};

    //-\Delat=f, g denotes the boundary condition
    EquationSolver(const int &_N, const Function &f){  
        N=_N;
        h=1.0/N;
        if(D==Domain::regular){
            for(int j=1; j<N; ++j){
                for(int i=1; i<N; ++i){
                    values.push_back(f(i*h, j*h)*h*h);
                    grids.push_back(vector<double>{i*h,j*h});
                }
            }
        }
        else{
            cerr<<"not regular domain"<<endl;
            return;
        }
    }
    EquationSolver(const int &_N, const Function &f, Circle *_c):N(_N), c(_c){  
        h=1.0/N;
        int count=0;
        for(int j=1; j<N; ++j){
            for(int i=1; i<N; ++i){
                grids.push_back(vector<double>{i*h,j*h});
                if(c->inCircle(i*h, j*h)){
                    values.push_back(0.0);
                    incircle.push_back(true);
                    count++;
                }
                else{
                    values.push_back(f(i*h, j*h)*h*h);
                    incircle.push_back(false);
                }
            }
        }
        if(count<4){
            cerr<<"invalid imput!"<<endl;
            return;
        }
    }

    vector<vector<double>> getgrids(){
        return grids;
    }

    vector<double> getvalues(){
        return values;
    }


    void norm_error(const Function &f,const string &filename){
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
        cout<<"------------------------------------------- calculating errors -----------------------------------------"<<endl;
        cout<<"l_1 norm is "<<l1_norm<<endl;
        cout<<"l_2 norm is "<<sqrt(l2_norm)<<endl;
        cout<<"l_infty norm is"<<infinity_norm<<endl;
        nlohmann::json j;
        j["boundary_condition"] = BC;  
        j["domain"] = D;              
        j["l1_norm"] = l1_norm;
        j["l2_norm"] = sqrt(l2_norm);
        j["linfinity_norm"] = infinity_norm;
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
        std::ofstream outFile(filename, std::ios::out | std::ios::trunc); 
        if (outFile.is_open()) {
            outFile << jsonDataArray.dump(4); 
            outFile.close();
        } else {
            cerr << "Error opening file " << filename << endl;
        }

    }

    void solveEquation(const Function &g, const double &initial=0.0, const vector<int> &mixed=vector<int>{0,0,0,0,0}){
        solve(g,mixed);
        this->adjust(initial);
    }

    void print(const string &filename, const Function &f){
        nlohmann::json j;
        j["boundary_condition"] = BC; 
        j["Domain"]=D;
        j["grids"] = grids;
        j["values_on_grids"] = values;
        vector<double> realvalues=this->realValues(f);
        j["real_values"]=values; 
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
        cout<<"data saved"<<endl;
    }
    
};






#endif