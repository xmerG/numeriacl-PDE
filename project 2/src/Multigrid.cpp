#include"Multigrid.h"
template<int dim>
Vector Multigrid<dim>::w_Jacobi(int j, const Vector &initial){
    Sparse_Matrix& A = discretors[j].first;
    double weighted=0.0;
    if constexpr(dim==1){
        weighted=1.0/3.0;
    }
    else {
        weighted=1.0/6.0;
    }
    Vector temp=A*initial;
    return initial+(discretors[j].second-temp)*weighted;
}

template<int dim>
Vector Multigrid<dim>::V_cycle(const int &_n, Vector& initial_guess, int nu1, int nu2){ //_n 网格个数
    for(int i=0; i<nu1; ++i){
        initial_guess=this->w_Jacobi(_n, initial_guess);
    }
    Vector v=initial_guess;
    if(_n==2){
        v=this->corsa_solve(_n);
    }
    else{
        Vector fnew=discretors[_n].second-discretors[_n].first*initial_guess;
        int corsa=_n/2;
        discretors[corsa].second=(*restriction)(fnew);
        if (dim==1){
            initial_guess.go_zero(corsa-1);
        }
        else{
            initial_guess.go_zero((corsa-1)*(corsa-1));
        }
        initial_guess=V_cycle(corsa, initial_guess, nu1, nu2);
        v=v+(*prolongation)(initial_guess);
    }
    for(int i=0; i<nu2; ++i){
        v=this->w_Jacobi(_n, v);
    }
    return v;    
}

template<int dim>
Vector Multigrid<dim>::FMG(const int &_n, int nu1, int nu2){
    int corsa=_n;
    if(_n==2){
        Vector v=Vector(1);
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

template<int dim>
Vector Multigrid<dim>::corsa_solve(const int &i){
    vector<double> matrix=discretors[i].first.convert_to_vector();
    vector<double> values=discretors[i].second.getelements();
    int n=discretors[i].first.getdim();
    vector<int> ipiv(n);
    int info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, 1, matrix.data(), n, ipiv.data(), values.data(), n);
    return Vector(values);
}

template<>
void Multigrid<1>::create_grids_D(const Function &f, const Function &g, const int &i) {
    Sparse_Matrix A(i-1);
    double h=1.0/i;
    Vector fh(i-1);
    for(int j=0; j<i-1; ++j){
        A.setValues(j,j,2.0);
        A.setValues(j,j-1,-1.0);
        A.setValues(j, j+1, -1.0);
    }
    if(i==n){
        for(int j=0; j<n-1; ++j){
            fh.set_Value(j, h*h*f((j+1)*h));
        }
        fh.set_Value(0, fh(0)+g(0.0));
        fh.set_Value(i-2, fh(i-2)+g(1.0));
    }
    discretors[i]=pair{move(A),move(fh)};   
}

template<>
void Multigrid<2>::create_grids_D(const Function &f, const Function &g, const int &i){
    int dim=pow(i-1,2);
    Sparse_Matrix A(dim);
    double h=1.0/i;
    Vector fh(dim);
    for(int j=0; j<dim; ++j){
        A.setValues(j,j, 4.0);
        A.setValues(j, j+1, -1.0);
        A.setValues(j, j-1, -1.0);
        A.setValues(j, j+i-1, -1.0);
        A.setValues(j, j-i+1, -1.0);
    }
    if(i==n){
        for(int j=0; j<i-1; ++j){
            for(int k=0; k<i-1; ++k){
                fh.set_Value(k, j, h*h*f((k+1)*h, (j+1)*h));
            }
        }
        for(int j=0; j<i-1; ++j){
            double temp=(j+1)*h;
            double value=fh(j, 0)+g(temp, 0.0);
            fh.set_Value(j,0, value);
            
            value=fh(0, j)+g(0.0, temp);
            fh.set_Value(0, j, value);
    
            value=fh(j, i-2)+g(temp, 1.0);
            fh.set_Value(j, i-2, value);
    
            value=fh(i-2, j)+g(1.0, temp);
            fh.set_Value(i-2, j, value);
        }
    }
    discretors[i]=make_pair(move(A), move(fh));
}




template<int dim>
Multigrid<dim>::Multigrid(){}

template<int dim>
Multigrid<dim>::Multigrid(const Function &f, const Function &g, BoundaryCondition bc,const int &i){
    n=i;
    if(bc==BoundaryCondition::Dirichlet){
        for(int j=i; j>1; j/=2){
            this->create_grids_D(f, g, j);
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