#include"Test.h"
#include"testFunction.h"
#include<fstream>
#include <chrono>
using json = nlohmann::json;

unique_ptr<Function> Laplacian(int dim){
    if(dim==1){
        return make_unique<Laplacian1>();
    }
    else{
        return make_unique<Laplacian2>();
    }
}

unique_ptr<Function> Primitive(int dim){
    if(dim==1){
        return make_unique<F1>();
    }
    else{
        return make_unique<F2>();
    }
}

unique_ptr<Function> BC_Function(int dim, string bc, vector<int> mix){
    if(dim==1 && bc=="Dirichlet"){
        return make_unique<F1>();
    }
    else if(dim==1 && bc=="Neumann"){
        return make_unique<Neumann1>();
    }
    else if(dim==1 &&bc=="Mixed"){
        return make_unique<Mixed1>(mix);  //-------------------------------------------
    }
    if(dim==2 && bc=="Dirichlet"){
        return make_unique<F2>();
    }
    else if(dim==2 && bc=="Neumann"){
        return make_unique<Neumann2>();
    }
    else{
        return make_unique<Mixed2>(mix);
    }

}

void write_test_results_to_csv(
    int n,
    const std::string& BC,
    const std::string& restriction_opr,
    const std::string& prolongation_opr,
    const std::string& Cycle,
    int dimension,
    const std::chrono::milliseconds& duration
) {
    std::ofstream csv_file("../output/test_results.csv", std::ios::app); 

    if (csv_file.is_open()) {
        if (csv_file.tellp() == 0) { 
            csv_file << "GridSize,BoundaryCondition,Restriction,Prolongationr,Cycle,Dimension,RunningTime(ms)\n";
        }

        csv_file << n << ","
                 << BC << ","
                 << restriction_opr << ","
                 << prolongation_opr << ","
                 << Cycle << ","
                 << dimension << ","
                 << duration.count() << "\n";

        csv_file.close();
        std::cout << "Test results saved to 'test_results.csv'" << std::endl;
    } else {
        std::cerr << "Failed to open CSV file for writing!" << std::endl;
    }
}

void test(const string &filename,  const string &output){
    ifstream inputFile(filename);
    json json_data;
    inputFile >> json_data;
    for(auto item : json_data){
        int n=item["grid_number"];
        BoundaryCondition bc;
        vector<int> mix=vector<int>{0,0,0,0};
        string BC=item["boundary_conditions"];
        if(BC=="Dirichlet"){
            bc=BoundaryCondition::Dirichlet;
        }
        else if(BC=="Neumann"){
            bc=BoundaryCondition::Neumann;
        }
        else{
            bc=BoundaryCondition::Mixed;
            if(item.contains("mixed")  && item["mixed"].is_array()){
                vector<int> v=item["mixed"];
                mix=v;
            }
        }
        string restriction_opr=item["restriction_operator"];
        string prolongation_opr=item["interpolation_operator"];
        string Cycle=item["cycle_type"];
        int max_itr_time=item["max_iterations"];
        double relative_accuracy=item["max_iterations"];
        double tol=item["relative_accuracy"];
        const int dimension=item["dimension"];
        int number=0;
        if(dimension==1){
            number=n+1;
        }
        else if(dimension==2){
            number=pow(n+1, 2);
        }
        Vector v=Vector(number);
        Vector &initial=v;
        unique_ptr<Function> laplacian_func = Laplacian(dimension);
        unique_ptr<Function> bc_func = BC_Function(dimension, BC, mix);
        unique_ptr<Function> f0=Primitive(dimension);
        double value=0.0;
        if(bc==BoundaryCondition::Neumann){
            if(dimension==1){
                value=(*f0)(0.0);
            }
            else{
                value=(*f0)(0.0, 0.0);
            }
        }
        auto start_time = chrono::high_resolution_clock::now();
        if(dimension==1){
            Multigrid<1> M(*laplacian_func, *bc_func, bc, n, mix);
            M.solve(restriction_opr, prolongation_opr, Cycle, initial, 5, 5, tol, value, max_itr_time);
            M.print_to_file(output, *f0);
        }
        else{
            Multigrid<2> M(*laplacian_func, *bc_func, bc, n, mix);
            M.solve(restriction_opr, prolongation_opr, Cycle, initial, 5, 5, tol, value, max_itr_time);
            M.print_to_file(output, *f0);
        }
        auto end_time = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
        cout<<"----------------------------------testing-------------------------------------------------------------------"<<endl;
        cout<<"grids: "<<n<<" "<<"Boundary condition: "<<BC<<" restrinction: "<<restriction_opr<<" prolongation: "<<prolongation_opr
            <<" cycle: "<<Cycle<<" dim: "<<dimension<<endl;
        cout<<"running time: "<<duration.count()<<"ms"<<endl;
        cout<<endl;
        write_test_results_to_csv(n, BC, restriction_opr, prolongation_opr, Cycle, dimension, duration);

    }
}