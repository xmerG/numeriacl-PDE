#include"Test.h"
#include"testFunction.h"
#include<fstream>
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



void Test(const string &filename){
    ifstream inputFile(filename);
    json json_data;
    inputFile >> json_data;
    for(auto item : json_data){
        int n=item["grid_number"];
        BoundaryCondition bc;
        vector<int> mix;
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
        int dimension=item["dimension"];
        int number=0;
        if(dimension==1){
            number=n+1;
        }
        else if(dimension==2){
            number=pow(n+1, 2);
        }
        Vector v=Vector(number);
        Vector &initial=v;

    }
}