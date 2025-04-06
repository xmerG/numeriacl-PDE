#include"Multigrid.h"
template<int dim>
Vector Multigrid<dim>::w_Jacobi(int j, const Vector &initial){
    Sparse_Matrix& A = discretors[j].first;
    double currh=1.0/j;
    double weighted=0.0;
    if constexpr(dim==1){
        weighted=1.0/3.0;
    }
    else {
        weighted=1.0/6.0;
    }
    Vector temp=A*initial;
    return initial+(discretors[j].second*(currh*currh)-temp)*weighted;
}


template<int dim>
Vector Multigrid<dim>::V_cycle(const int &_n, Vector& initial_guess, int nu1, int nu2){ //_n 网格个数
    for(int i=0; i<nu1; ++i){
        initial_guess=this->w_Jacobi(_n, initial_guess);
    }
    if(_n==4){
        Vector v=discretors[_n].second*(1.0/(_n*_n));
        for(int i=0; i<10; ++i){
            discretors[_n].first.Gauss_Seidel(initial_guess, v);
        }
        for(int i=0; i<nu2; ++i){
            initial_guess=this->w_Jacobi(_n, initial_guess);
        }
        return initial_guess;
    }
    else{
        Vector fnew=discretors[_n].second-discretors[_n].first*initial_guess*(_n*_n);
        int corsa=_n/2;
        discretors[corsa].second=(*restriction)(fnew);
        Vector v(corsa-1);
        if (dim==2){
            v.go_zero((corsa-1)*(corsa-1));
        }
        v=V_cycle(corsa,v , nu1, nu2);
        initial_guess=initial_guess+(*prolongation)(v);
        for(int i=0; i<nu2; ++i){
            initial_guess=this->w_Jacobi(_n, initial_guess);
            //discretors[_n].first.Gauss_Seidel(initial_guess, discretors[_n].second);
        }
        return initial_guess;    
    }
}

template<int dim>
Vector Multigrid<dim>::FMG(const int &_n, int nu1, int nu2){
    int corsa=_n;
    if(_n==4){
        Vector v=Vector(_n-1);
        v=this->V_cycle(_n, v, nu1, nu2);
        return v;
    }
    else{
        corsa/=2;
        discretors[corsa].second=(*restriction)(discretors[_n].second);
        Vector v=this->FMG(corsa, nu1, nu2);
        v=(*prolongation)(v);
        v=this->V_cycle(_n, v, nu1, nu2);
        return v;
    }
}

/*template<int dim>
vector<double> Multigrid<dim>::corsa_solve(const int &i) {
    vector<double> matrix = discretors[i].first.convert_to_vector();
    vector<double> values = discretors[i].second.getelements();
    int n = discretors[i].first.getdim();
    vector<int> ipiv(n);
    int info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, 1, matrix.data(), n, ipiv.data(), values.data(), n);
    return values;
}*/

template<>
void Multigrid<1>::create_grids_D(const Function &f, const Function &g, const int &i) {
    Sparse_Matrix A(i+1);
    Vector fh(i+1);
    for(int j=1; j<i; ++j){
        A.setValues(j,j,2.0);
        A.setValues(j,j-1,-1.0);
        A.setValues(j, j+1, -1.0);
    }
    A.setValues(0,0, 2.0);
    A.setValues(i, i, 2.0);
    if(i==n){
        double h=1.0/i;
        for(int j=1; j<n; ++j){
            fh.set_Value(j,f(j*h));
        }
        fh.set_Value(0, 2*g(0.0)*n*n);
        fh.set_Value(n, 2*g(1.0)*n*n);
    }
    discretors[i]=pair{move(A),move(fh)};   
}

template<>
void Multigrid<2>::create_grids_D(const Function &f, const Function &g, const int &i){
    int dim=pow(i+1,2);
    Sparse_Matrix A(dim);
    double h=1.0/i;
    Vector fh(dim);
    //内部
    for(int j=1; j<i; ++j){
        for(int k=1; k<i; ++k){
            int index=k+j*(i+1);
            A.setValues(index, index, 4.0);
            A.setValues(index, index+1, -1.0);
            A.setValues(index, index-1, -1.0);
            A.setValues(index, index+i-1, -1.0);
            A.setValues(index, index-i+1, -1.0);
        }
    }

    //设置边界值
    for(int j=0; j<=dim; ++j){
        A.setValues(j, j, 4.0);
        int index=j*(i+1);
        A.setValues(index, index, 4.0);
        index+=i;
        A.setValues(index, index, 4.0);
        index=(i+1)*i+j;
        A.setValues(index, index, 4.0);
        /*A.setValues(j,j, 4.0);
        A.setValues(j, j+1, -1.0);
        A.setValues(j, j-1, -1.0);
        A.setValues(j, j+i-1, -1.0);
        int index=j*(i-1);
        A.setValues(index, index, 4.0);
        A.setValues(index, index+1, -1.0);
        A.setValues(index, index+i-1, -1.0);
        A.setValues(index, index-i+1, -1.0);
        index=index+i-2;
        A.setValues(index, index, 4.0);
        A.setValues(index, index-1, -1.0);
        A.setValues(index, index-i+1, -1.0);
        A.setValues(index, index+i-1, -1.0);
        index=(i-1)*(i-2)+j;
        A.setValues(index, index, 4.0);
        A.setValues(index, index-1, -1.0);
        A.setValues(index, index+1, -1.0);
        A.setValues(index, index-i+1, -1.0);*/
    }
    /*A.setValues(0,0, 4.0);
    A.setValues(0,1, -1.0);
    A.setValues(0, i-1, -1.0);
    A.setValues(i-2, i-2, 4.0);
    A.setValues(i-2, i-3, -1.0);
    A.setValues(i-2, 2*i-3, -1.0);
    A.setValues(dim-1, dim-1, 4.0);
    A.setValues(dim-1, dim-2, -1.0);
    A.setValues(dim-1, dim-i+1, -1.0);
    int index=(i-1)*(i-2);
    A.setValues(index, index, 4.0);
    A.setValues(index, index+1, -1.0);
    A.setValues(index, index-i+1, -1.0);*/
    if(i==n){
        int n2=n*n;
        for(int j=1; j<n; ++j){
            for(int k=1; k<n; ++k){
                fh.set_Value(k, j, f(k*h, j*h));
            } 
        }
        //边界
        for(int j=0; j<=n; ++j){
            double temp=j*h;
            fh.set_Value(j,0, g(temp, 0.0)*n2);
            
            fh.set_Value(0, j, g(0.0, temp)*n2);
    
            fh.set_Value(j, n, g(temp, 1.0)*n2);
    
            fh.set_Value(n, j, g(1.0, temp)*n2);
        }
    }
    discretors[i]=make_pair(move(A), move(fh));
}

//---------------------------------------to be modified--------------------------------------------------------------
template<>
void Multigrid<1>::create_grids_N(const Function &f, const Function &g,const int &i, double value){
    Sparse_Matrix A(i-1);
    Vector fh(i-1);
    for(int j=1; j<i-2; ++j){
        A.setValues(j,j,2.0);
        A.setValues(j,j-1,-1.0);
        A.setValues(j, j+1, -1.0);
    }
    A.setValues(0,0, 2.0);
    A.setValues(i-2, i-3, -2.0);
    A.setValues(i-2, i-2, 2.0);
    if(i==n){
        double h=1.0/i;
        for(int j=0; j<n-1; ++j){
            fh.set_Value(j,f((j+1)*h));
        }
        fh.set_Value(0, 2*value*n*n);
        fh.set_Value(i-2, 3*fh(i-2)+2*g(1.0)*n);
    }
    discretors[i]=pair{move(A),move(fh)};   
}

//-------------------------------------------------to be modified----------------------------------------
template<>
void Multigrid<2>::create_grids_N(const Function &f, const Function &g,const int &i, double value){
    int dim=pow(i-1,2);
    Sparse_Matrix A(dim);
    double h=1.0/i;
    Vector fh(dim);
    int temp=(i-1)*(i-2);
    for(int j=1; j<dim; ++j){
        A.setValues(j,j, 4.0);
        if(j==i-2 || j==dim-1 || i==temp){
            A.setValues(j, j+1, -2.0);
            A.setValues(j, j-1, -2.0);
            A.setValues(j, j+i-1, -2.0);
            A.setValues(j, j-i+1, -2.0);
        }
        else{
            A.setValues(j, j+1, -1.0);
            A.setValues(j, j-1, -1.0);
            A.setValues(j, j+i-1, -1.0);
            A.setValues(j, j-i+1, -1.0);
        }
    }
    for(int j=1; j<i-2; ++j){
        A.setValues(j, j+i-1, -1.0);
        A.setValues(j, j-1, -1.5);
        A.setValues(j, j+1, -1.5);
        int index=j*(i-1);
        A.setValues(index, index+1, -1.0);
        A.setValues(index, index+i-1, -1.5);
        A.setValues(index, index-i+1, -1.5);
        index+=i-2;
        A.setValues(index, index-1, -1.0);
        A.setValues(index, index+i-1, -1.5);
        A.setValues(index, index-i+1, -1.5);
        index=(i-1)*(i-2)+j;
        A.setValues(index, index+i+1,-1.0);
        A.setValues(index, index-1, -1.5);
        A.setValues(index, index+1, -1.5);
        A.setValues(0, 0, 4.0);
    }
    if(i==n){
        for(int j=0; j<i-1; ++j){
            for(int k=0; k<i-1; ++k){
                fh.set_Value(k, j, f((k+1)*h, (j+1)*h));
            }
        }
        for(int j=1; j<i-2; ++j){
            double temp=(j+1)*h;
            double value=1.5*fh(j, 0)+g(temp, 0.0)*n;
            fh.set_Value(j,0, value);
            
            value=1.5*fh(0, j)+g(0.0, temp)*n;
            fh.set_Value(0, j, value);
    
            value=1.5*fh(j, i-2)+g(temp, 1.0)*n;
            fh.set_Value(j, i-2, value);
    
            value=1.5*fh(i-2, j)+g(1.0, temp)*n;
            fh.set_Value(i-2, j, value);
        }
        fh.set_Value(0,0, 4*value*n*n);
        fh.set_Value(i-2,0, 3*fh(i-2,0)+(g(1-h,0.0)+g(1.0, h))*n*2);
        fh.set_Value(0,i-2, 3*fh(0,i-2)+(g(h,1.0)+g(0.0, 1-h))*n*2);
        fh.set_Value(i-2,i-2, 3*fh(i-2,i-2)+(g(1-h,1.0)+g(1.0, 1-h))*n*2);
    }
    discretors[i]=make_pair(move(A), move(fh));
}


template<int dim>
Multigrid<dim>::Multigrid(){}

template<int dim>
Multigrid<dim>::Multigrid(const Function &f, const Function &g, BoundaryCondition bc,const int &i, double value){
    n=i;
    if(bc==BoundaryCondition::Dirichlet){
        for(int j=i; j>3; j/=2){
            this->create_grids_D(f, g, j);
        }
    }
    else if(bc==BoundaryCondition::Neumann){
        for(int j=i; j>3; j/=2){
            this->create_grids_N(f, g, j, value);
        }
    }
}

template<int dim>
void Multigrid<dim>::print(){
    for(auto& entry : discretors){
        int j=entry.first;
        cout<<j<<endl;
        couplet& c = entry.second;
        c.first.print();
        cout<<endl;
    }
}

template<int dim>
void Multigrid<dim>::solve(const string& r, const string& p, const string& c, 
    Vector& initial_guess, int nu1, int nu2) {
    if (r == "full_weighting") {
    restriction = make_unique<Full_weighting<dim>>(); 
    } 
    else if (r == "injection") {
    restriction = make_unique<Injection<dim>>(); 
    }
     else {
    cerr << "Invalid restriction method: " << r << endl;
    return;
    }
    
    if (p == "linear") {
    prolongation = make_unique<Linear<dim>>(); 
    } 
    else if (p == "quadric") { 
    prolongation = make_unique<Quadric<dim>>();  
    }
     else {
    cerr << "Invalid prolongation method: " << p << endl;
    return;
    }
    if(c=="v-cylce"){
        solutions=this->V_cycle(n, initial_guess, nu1, nu2);
        Vector f=discretors[n].second;
        int n2=n*n;
        for(int i=0; i<5; ++i){
            Vector error=f-discretors[n].first*solutions*n2;
            discretors[n].second=move(error);
            solutions=solutions+this->V_cycle(n, initial_guess, nu1, nu2);
        }
    }
    else{
        solutions=this->FMG(n, nu1, nu2);
    }
}

template<>
vector<double> Multigrid<1>::error(const Function &f){
    double h=1.0/n;
    vector<double> e(n-1,0.0);
    for(int j=0; j<n-1; ++j){
        e[j]=abs(solutions(j)-f((j+1)*h));
    }
    return e;
}

template<>
vector<double> Multigrid<2>::error(const Function &f){
    double h=1.0/n;
    int m=(n-1)*(n-1);
    vector<double> e(m, 0.0);
    for(int i=0; i<n-1; ++i){
        for(int j=0; j<n-1; ++j){
            e[i*(n-1)+j]=abs(solutions(j, i)-f((j+1)*h, (i+1)*h));
        }
    }
    return e;
}



template<int dim>
void Multigrid<dim>::print_to_file(const string &filename, BoundaryCondition BC,const Function &f) {
    nlohmann::json j;
    j["boundary_condition"] = BC; 
    vector<double> e=this->error(f);//-------------------------------------------------------;
    j["errors"] = e;
    j["solutions"]=solutions.getelements();
    j["number"]=n;
    ifstream file_check(filename); 
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


template class Multigrid<1>;
template class Multigrid<2>;